!##############################################################################
!# ****************************************************************************
!# <name> element_quad2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the implementations of the 2D quadrilateral basis
!# functions.
!#
!# </purpose>
!##############################################################################

module element_quad2d

!$use omp_lib
  use basicgeometry
  use derivatives
  use elementbase
  use fsystem
  use mprimitives
  use perfconfig
  use transformation

  implicit none

  private

  public :: elem_Q0
  public :: elem_Q0_mult
  public :: elem_Q0_sim
  public :: elem_Q1
  public :: elem_Q1_mult
  public :: elem_Q1_sim
  public :: elem_Q2
  public :: elem_Q2_mult
  public :: elem_Q2_sim
  public :: elem_QP1
  public :: elem_QP1NP
  public :: elem_EM30
  public :: elem_EM30_mult
  public :: elem_EM30_sim
  public :: elem_E030
  public :: elem_E030_mult
  public :: elem_E030_sim
  public :: elem_EB30
  public :: elem_EB30_mult
  public :: elem_EB30_sim
  public :: elem_EM31
  public :: elem_EM31_mult
  public :: elem_EM31_sim
  public :: elem_E031
  public :: elem_E031_mult
  public :: elem_E031_sim
  public :: elem_E050
  public :: elem_E050_mult
  public :: elem_E050_sim
  public :: elem_EB50
  public :: elem_EB50_mult
  public :: elem_EB50_sim
  public :: elem_eval_Q1_2D
  public :: elem_eval_Q1B_2D
  public :: elem_eval_EM11_2D
  public :: elem_eval_Q2_2D
  public :: elem_eval_Q2H_2D
  public :: elem_eval_Q3_2D
  public :: elem_eval_QP1_2D
  public :: elem_eval_QPW4P0_2D
  public :: elem_eval_QPW4P1_2D
  public :: elem_eval_QPW4P1T_2D
  public :: elem_eval_QPW4P2_2D
  public :: elem_eval_QPW4DCP1_2D
  public :: elem_eval_QPW4P1TVDF_2D
  public :: elem_eval_E030_2D
  public :: elem_eval_EB30_2D
  public :: elem_eval_Q1TBNP_2D
  public :: elem_eval_EN30_2D
  public :: elem_eval_EN31_2D
  public :: elem_eval_E032_2D
  public :: elem_eval_E050_2D
  public :: elem_eval_EB50_2D
  public :: elem_eval_EN50_2D
  public :: elem_eval_EN51_2D
  public :: elem_eval_DCQP1_2D
  public :: elem_eval_DCQP2_2D
  public :: elem_DG_T0_2D
  public :: elem_DG_T0_2D_mult
  public :: elem_DG_T0_2D_sim
  public :: elem_DG_T1_2D
  public :: elem_DG_T1_2D_mult
  public :: elem_DG_T1_2D_sim
  public :: elem_DG_T2_2D
  public :: elem_DG_T2_2D_mult
  public :: elem_DG_T2_2D_sim
  public :: elem_DG_T3_2D
  public :: elem_DG_T3_2D_mult
  public :: elem_DG_T3_2D_sim

contains

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

  pure subroutine elem_Q0 (celement, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q0.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_Q0_mult (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0.
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
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints

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

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_Q0_sim (celement, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
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
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  !</input>

  !<output>

  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th
  !   basis function of the finite element in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas

  !</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  end subroutine

!**************************************************************************
! Element subroutines for parametric Q1 element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_Q1 (celement, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D) :: Dhelp
  real(DP) :: dx,dy

  real(DP) :: dxj !auxiliary variable

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
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  ! If function values are desired, calculate them.
  !
  ! We have to compute the basis functions in the points.
  ! I.e., we have to evaluate
  !
  !    phi_k(x) = Pk(sigma^-1(x)) = Pk^(x^)
  !
  ! with x being the real world coordinates, x^ the coordinates
  ! in the reference element and sigma: x^ -> x the mapping
  ! between the reference and the real element.
  ! So just evaluate the basis functions in the points x^ 
  ! on the reference element

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
  !
  ! We have to evaluate "grad(phi(x))". This is done by
  ! using the chain rule as follows:
  !
  !    grad(phi_k(x))^T = D phi_k(x)
  !                     = D Pk(sigma^-1(x))
  !                     = D Pk(x^) * D sigma^-1(x)
  !
  ! Now note that the Jacobian "D sigma(x^)" of the mapping
  ! between the reference and the real element is given in Djac:
  !
  !    D sigma(x^) = ( Djac(1) Djac(3) )
  !                  ( Djac(2) Djac(4) )
  !
  ! Its inverse then reads
  !
  !    D sigma^-1(x) = 1/det ( -Djac(4)  Djac(3) )
  !                          (  Djac(2) -Djac(1) )
  !
  ! with det = Djac(1)*Djac(4) - Djac(2)*Djac(3).
  ! So all in all, the derivative can be computed as follows:
  !
  !    grad(phi_k(x))^T = D Pk(x^) * D sigma^-1(x)
  !                     = ( dx(Pk(x^)) dy(Pk(x^)) ) * 1/det ( -Djac(4)  Djac(3) )
  !                                                         (  Djac(2) -Djac(1) )
  ! i.e., we have
  !
  !   dx(phi_k(x)) = 1/det ( -dx(Pk(x^))*Djac(4) + dy(Pk(x^))*Djac(2) )
  !   dy(phi_k(x)) = 1/det (  dx(Pk(x^))*Djac(3) - dy(Pk(x^))*Djac(1) )
  !

!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    ! We take ddet = 0.25/det and thus, put the constant factor 0.25
    ! in front of the basis functions into this multiplier.
    ! Dhelp receives the derivatives dx(Pk(x^)) and dy(Pk(x^))
    ! without the constant factor 0.25 (saves some multiplications).

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

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_Q1_mult (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC,i) = 0.25E0_DP*(1E0_DP-Dpoints(1,i))*(1E0_DP-Dpoints(2,i))
      Dbas(2,DER_FUNC,i) = 0.25E0_DP*(1E0_DP+Dpoints(1,i))*(1E0_DP-Dpoints(2,i))
      Dbas(3,DER_FUNC,i) = 0.25E0_DP*(1E0_DP+Dpoints(1,i))*(1E0_DP+Dpoints(2,i))
      Dbas(4,DER_FUNC,i) = 0.25E0_DP*(1E0_DP-Dpoints(1,i))*(1E0_DP+Dpoints(2,i))
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj(:) = 0.25E0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
      Dhelp(1,1,i) =-(1E0_DP-Dpoints(2,i))
      Dhelp(2,1,i) = (1E0_DP-Dpoints(2,i))
      Dhelp(3,1,i) = (1E0_DP+Dpoints(2,i))
      Dhelp(4,1,i) =-(1E0_DP+Dpoints(2,i))
      Dhelp(1,2,i) =-(1E0_DP-Dpoints(1,i))
      Dhelp(2,2,i) =-(1E0_DP+Dpoints(1,i))
      Dhelp(3,2,i) = (1E0_DP+Dpoints(1,i))
      Dhelp(4,2,i) = (1E0_DP-Dpoints(1,i))
    end do

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_Q1_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, &
                          Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 0.25_DP*(1.0_DP-Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dbas(2,DER_FUNC,i,j) = 0.25_DP*(1.0_DP+Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dbas(3,DER_FUNC,i,j) = 0.25_DP*(1.0_DP+Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
        Dbas(4,DER_FUNC,i,j) = 0.25_DP*(1.0_DP-Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,Dxj,Dhelp) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 0.25E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
        Dhelp(1,1,i) =-(1E0_DP-Dpoints(2,i,j))
        Dhelp(2,1,i) = (1E0_DP-Dpoints(2,i,j))
        Dhelp(3,1,i) = (1E0_DP+Dpoints(2,i,j))
        Dhelp(4,1,i) =-(1E0_DP+Dpoints(2,i,j))
        Dhelp(1,2,i) =-(1E0_DP-Dpoints(1,i,j))
        Dhelp(2,2,i) =-(1E0_DP+Dpoints(1,i,j))
        Dhelp(3,2,i) = (1E0_DP+Dpoints(1,i,j))
        Dhelp(4,2,i) = (1E0_DP-Dpoints(1,i,j))
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(1,1,i) &
                                    - Djac(2,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(2,1,i) &
                                    - Djac(2,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(3,1,i) &
                                    - Djac(2,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(4,1,i) &
                                    - Djac(2,i,j) * Dhelp(4,2,i))
        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(1,1,i) &
                                    - Djac(1,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(2,1,i) &
                                    - Djac(1,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(3,1,i) &
                                    - Djac(1,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(4,1,i) &
                                    - Djac(1,i,j) * Dhelp(4,2,i))
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric Q2 element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_Q2 (celement, Dcoords, Djac, ddetj, Bder, &
                      Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(9,5) :: Dhelp
  real(DP) :: dx,dy

  integer :: idof
  real(DP) :: dxj,dxjs,dc1,dc2,dc3 !auxiliary variable

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q4 = 0.25_DP

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
      do idof = 1,9
        Dbas(idof,DER_DERIV_X) = &
            dxj * (Djac(4) * Dhelp(idof,1) - Djac(2) * Dhelp(idof,2))
      end do
    endif

    ! y-derivatives on current element
    if (Bder(DER_DERIV_Y)) then
      do idof = 1,9
        Dbas(idof,DER_DERIV_Y) = &
            -dxj * (Djac(3) * Dhelp(idof,1) - Djac(1) * Dhelp(idof,2))
      end do
    endif
  endif

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then
    dxj = 1.0_DP / ddetj

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]^2->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]
    !

    dxjs = dxj*dxj

    !xx-, xy and yy-derivatives on reference element
    Dhelp(1,3) = Q2*(-1.0_DP+dy)*dy
    Dhelp(2,3) = Q2*(-1.0_DP+dy)*dy
    Dhelp(3,3) = Q2*(1.0_DP+dy)*dy
    Dhelp(4,3) = Q2*(1.0_DP+dy)*dy
    Dhelp(5,3) = -(-1.0_DP+dy)*dy
    Dhelp(6,3) = 1.0_DP-dy**2
    Dhelp(7,3) = -(1.0_DP+dy)*dy
    Dhelp(8,3) = 1.0_DP-dy**2
    Dhelp(9,3) = -2.0_DP+2.0_DP*dy**2

    Dhelp(1,4) = dx*dy-Q2*dx-Q2*dy+Q4
    Dhelp(2,4) = dx*dy-Q2*dx+Q2*dy-Q4
    Dhelp(3,4) = dx*dy+Q2*dx+Q2*dy+Q4
    Dhelp(4,4) = dx*dy+Q2*dx-Q2*dy-Q4
    Dhelp(5,4) = -2.0_DP*dx*dy+dx
    Dhelp(6,4) = -2.0_DP*dx*dy-dy
    Dhelp(7,4) = -2.0_DP*dx*dy-dx
    Dhelp(8,4) = -2.0_DP*dx*dy+dy
    Dhelp(9,4) = 4.0_DP*dx*dy

    Dhelp(1,5) = Q2*(-1.0_DP+dx)*dx
    Dhelp(2,5) = Q2*(1.0_DP+dx)*dx
    Dhelp(3,5) = Q2*(1.0_DP+dx)*dx
    Dhelp(4,5) = Q2*(-1.0_DP+dx)*dx
    Dhelp(5,5) = 1.0_DP-dx**2
    Dhelp(6,5) = -(1.0_DP+dx)*dx
    Dhelp(7,5) = 1.0_DP-dx**2
    Dhelp(8,5) = -(-1.0_DP+dx)*dx
    Dhelp(9,5) = -2.0_DP+2.0_DP*dx**2

    ! WARNING: NOT TESTED!!!

    ! xx-derivatives on current element
    if (Bder(DER_DERIV2D_XX)) then
      dc1 = dxjs * Djac(4)**2
      dc2 = dxjs * ( -2.0_DP * Djac(4) * Djac(2) )
      dc3 = dxjs * Djac(2)**2
      do idof = 1,9
        Dbas(idof,DER_DERIV2D_XX) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif

    ! xy-derivatives on current element
    if (Bder(DER_DERIV2D_XY)) then
      dc1 = - dxjs * Djac(4) * Djac(3)
      dc2 = dxjs * ( Djac(4) * Djac(1) + Djac(3) * Djac(2) )
      dc3 = - dxjs * Djac(1) * Djac(4)
      do idof = 1,9
        Dbas(idof,DER_DERIV2D_XY) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif

    ! yy-derivatives on current element
    if (Bder(DER_DERIV2D_YY)) then
      dc1 = dxjs * Djac(3)**2
      dc2 = dxjs * ( -2.0_DP * Djac(3) * Djac(1) )
      dc3 = dxjs * Djac(1)**2
      do idof = 1,9
        Dbas(idof,DER_DERIV2D_YY) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif
  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_Q2_mult (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(9,5,npoints) :: Dhelp
  real(DP) :: dx,dy
  real(DP),dimension(npoints) :: Dxj,Dxjs,Dc1,Dc2,Dc3 !auxiliary variable

  integer :: i,idof   ! point counter

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q4 = 0.25_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
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
    end do
  !ENDIF

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    Dxj(:) = 1.0E0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
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
    end do

    ! x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        do idof = 1,9
          Dbas(idof,DER_DERIV_X,i) = &
              Dxj(i) * (Djac(4,i) * Dhelp(idof,1,i) - Djac(2,i) * Dhelp(idof,2,i))
        end do
!      END DO
!    ENDIF

    ! y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        do idof = 1,9
          Dbas(idof,DER_DERIV_Y,i) = &
              -Dxj(i) * (Djac(3,i) * Dhelp(idof,1,i) - Djac(1,i) * Dhelp(idof,2,i))
        end do
      end do
!    ENDIF
  endif

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]

    !x- and y-derivatives on reference element
    Dxj(:) = 1.0E0_DP / Ddetj(1:npoints)
    Dxjs = Dxj*Dxj

    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)

      !xx-, xy and yy-derivatives on reference element
      Dhelp(1,3,i) = Q2*(-1.0_DP+dy)*dy
      Dhelp(2,3,i) = Q2*(-1.0_DP+dy)*dy
      Dhelp(3,3,i) = Q2*(1.0_DP+dy)*dy
      Dhelp(4,3,i) = Q2*(1.0_DP+dy)*dy
      Dhelp(5,3,i) = -(-1.0_DP+dy)*dy
      Dhelp(6,3,i) = 1.0_DP-dy**2
      Dhelp(7,3,i) = -(1.0_DP+dy)*dy
      Dhelp(8,3,i) = 1.0_DP-dy**2
      Dhelp(9,3,i) = -2.0_DP+2.0_DP*dy**2

      Dhelp(1,4,i) = dx*dy-Q2*dx-Q2*dy+Q4
      Dhelp(2,4,i) = dx*dy-Q2*dx+Q2*dy-Q4
      Dhelp(3,4,i) = dx*dy+Q2*dx+Q2*dy+Q4
      Dhelp(4,4,i) = dx*dy+Q2*dx-Q2*dy-Q4
      Dhelp(5,4,i) = -2.0_DP*dx*dy+dx
      Dhelp(6,4,i) = -2.0_DP*dx*dy-dy
      Dhelp(7,4,i) = -2.0_DP*dx*dy-dx
      Dhelp(8,4,i) = -2.0_DP*dx*dy+dy
      Dhelp(9,4,i) = 4.0_DP*dx*dy

      Dhelp(1,5,i) = Q2*(-1.0_DP+dx)*dx
      Dhelp(2,5,i) = Q2*(1.0_DP+dx)*dx
      Dhelp(3,5,i) = Q2*(1.0_DP+dx)*dx
      Dhelp(4,5,i) = Q2*(-1.0_DP+dx)*dx
      Dhelp(5,5,i) = 1.0_DP-dx**2
      Dhelp(6,5,i) = -(1.0_DP+dx)*dx
      Dhelp(7,5,i) = 1.0_DP-dx**2
      Dhelp(8,5,i) = -(-1.0_DP+dx)*dx
      Dhelp(9,5,i) = -2.0_DP+2.0_DP*dx**2

    end do

    ! WARNING: NOT TESTED!!!

    ! xx-derivatives on current element
    !if (Bder(DER_DERIV_XX)) then
      do i=1,npoints
        Dc1(i) = Dxjs(i) * Djac(4,i)**2
        Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(4,i) * Djac(2,i) )
        Dc3(i) = Dxjs(i) * Djac(2,i)**2
      end do

      do i=1,npoints
        do idof = 1,9
          Dbas(idof,DER_DERIV_XX,i) = &
              Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
        end do
      end do
    !endif

    ! xy-derivatives on current element
    !if (Bder(DER_DERIV_XY)) then
      do i=1,npoints
        Dc1(i) = - Dxjs(i) * Djac(4,i) * Djac(3,i)
        Dc2(i) = Dxjs(i) * ( Djac(4,i) * Djac(1,i) + Djac(3,i) * Djac(2,i) )
        Dc3(i) = - Dxjs(i) * Djac(1,i) * Djac(4,i)
      end do

      do i=1,npoints
        do idof = 1,9
          Dbas(idof,DER_DERIV_XY,i) = &
              Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
        end do
      end do
    !endif

    ! yy-derivatives on current element
    !if (Bder(DER_DERIV_YY)) then
      do i=1,npoints
        Dc1(i) = Dxjs(i) * Djac(3,i)**2
        Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(3,i) * Djac(1,i) )
        Dc3(i) = Dxjs(i) * Djac(1,i)**2
      end do

      do idof = 1,9
        Dbas(idof,DER_DERIV_YY,i) = &
            Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
      end do
    !endif
  endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_Q2_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, &
                          Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(9,5,npoints) :: Dhelp
  real(DP) :: dx,dy
  real(DP),dimension(npoints) :: Dxj,Dxjs,Dc1,Dc2,Dc3 !auxiliary variables

  integer :: i   ! point counter
  integer :: j   ! element counter
  integer :: idof

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q4 = 0.25_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i,dx,dy) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
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
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,Dxj,dx,dy,Dhelp,idof) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1.0E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
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
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          do idof = 1,9
            Dbas(idof,DER_DERIV_X,i,j) = &
                Dxj(i) * (Djac(4,i,j) * Dhelp(idof,1,i) &
                          - Djac(2,i,j) * Dhelp(idof,2,i))
!          end do
!        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
!        do i=1,npoints
!          do idof = 1,9
            Dbas(idof,DER_DERIV_Y,i,j) = &
                -Dxj(i) * (Djac(3,i,j) * Dhelp(idof,1,i) &
                - Djac(1,i,j) * Dhelp(idof,2,i))
          end do
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]

    !x- and y-derivatives on reference element
    !$omp parallel do default(shared)&
    !$omp private(i,dx,dy,Dxj,Dxjs,Dhelp,Dc1,Dc2,Dc3,idof)&
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1.0E0_DP / Ddetj(1:npoints,j)
      Dxjs = Dxj*Dxj

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)

        !xx-, xy and yy-derivatives on reference element
        Dhelp(1,3,i) = Q2*(-1.0_DP+dy)*dy
        Dhelp(2,3,i) = Q2*(-1.0_DP+dy)*dy
        Dhelp(3,3,i) = Q2*(1.0_DP+dy)*dy
        Dhelp(4,3,i) = Q2*(1.0_DP+dy)*dy
        Dhelp(5,3,i) = -(-1.0_DP+dy)*dy
        Dhelp(6,3,i) = 1.0_DP-dy**2
        Dhelp(7,3,i) = -(1.0_DP+dy)*dy
        Dhelp(8,3,i) = 1.0_DP-dy**2
        Dhelp(9,3,i) = -2.0_DP+2.0_DP*dy**2

        Dhelp(1,4,i) = dx*dy-Q2*dx-Q2*dy+Q4
        Dhelp(2,4,i) = dx*dy-Q2*dx+Q2*dy-Q4
        Dhelp(3,4,i) = dx*dy+Q2*dx+Q2*dy+Q4
        Dhelp(4,4,i) = dx*dy+Q2*dx-Q2*dy-Q4
        Dhelp(5,4,i) = -2.0_DP*dx*dy+dx
        Dhelp(6,4,i) = -2.0_DP*dx*dy-dy
        Dhelp(7,4,i) = -2.0_DP*dx*dy-dx
        Dhelp(8,4,i) = -2.0_DP*dx*dy+dy
        Dhelp(9,4,i) = 4.0_DP*dx*dy

        Dhelp(1,5,i) = Q2*(-1.0_DP+dx)*dx
        Dhelp(2,5,i) = Q2*(1.0_DP+dx)*dx
        Dhelp(3,5,i) = Q2*(1.0_DP+dx)*dx
        Dhelp(4,5,i) = Q2*(-1.0_DP+dx)*dx
        Dhelp(5,5,i) = 1.0_DP-dx**2
        Dhelp(6,5,i) = -(1.0_DP+dx)*dx
        Dhelp(7,5,i) = 1.0_DP-dx**2
        Dhelp(8,5,i) = -(-1.0_DP+dx)*dx
        Dhelp(9,5,i) = -2.0_DP+2.0_DP*dx**2

      end do

      ! WARNING: NOT TESTED!!!

      ! xx-derivatives on current element
      !if (Bder(DER_DERIV_XX)) then
        do i=1,npoints
          Dc1(i) = Dxjs(i) * Djac(4,i,j)**2
          Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(4,i,j) * Djac(2,i,j) )
          Dc3(i) = Dxjs(i) * Djac(2,i,j)**2
        end do

        do i=1,npoints
          do idof = 1,9
            Dbas(idof,DER_DERIV_XX,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

      ! xy-derivatives on current element
      !if (Bder(DER_DERIV_XY)) then
        do i=1,npoints
          Dc1(i) = - Dxjs(i) * Djac(4,i,j) * Djac(3,i,j)
          Dc2(i) = Dxjs(i) * ( Djac(4,i,j) * Djac(1,i,j) + Djac(3,i,j) * Djac(2,i,j) )
          Dc3(i) = - Dxjs(i) * Djac(1,i,j) * Djac(4,i,j)
        end do

        do i=1,npoints
          do idof = 1,9
            Dbas(idof,DER_DERIV_XY,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

      ! yy-derivatives on current element
      !if (Bder(DER_DERIV_YY)) then
        do i=1,npoints
          Dc1(i) = Dxjs(i) * Djac(3,i,j)**2
          Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(3,i,j) * Djac(1,i,j) )
          Dc3(i) = Dxjs(i) * Djac(1,i,j)**2
        end do

        do i=1,npoints
          do idof = 1,9
            Dbas(idof,DER_DERIV_YY,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric QP1 element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_QP1 (celement, Dcoords, Djac, ddetj, Bder, &
                            Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_QP1.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! local variables
  ! REAL(DP), DIMENSION(2,2) :: Dhelp
  real(DP) :: dxj

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

    dxj = 1.0_DP/ddetj

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

  end subroutine

!**************************************************************************
! Element subroutines for nonparametric QP1 element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_QP1NP (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_QP1NP.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

    ! The element is given the function value and the X- and Y-derivative
    ! in the midpoint of the reference element.
    ! That means: p(z1,z2) = a*1 + b*z1 + c*z2, coefficients (a,b,c).
    !
    ! The element is nonparametric, i.e. the coordinate system is element
    ! dependent. The point (z1,z2) = K(x1,x2) is the representation of the
    ! point (x1,x2) in the new coordinate system. The mapping K() maps
    ! a (x1,x2) to the normalised coordinate system defined by the midpoints
    ! of the edges.

    ! This element clearly works only with standard quadrilaterals
    integer, parameter :: NVE = 4

    ! auxiliary variables
    real(DP) :: dx, dy, dxj, dmx, dmy
    real(DP),dimension(4) :: DXM,DYM
    real(dp) :: deta1,dxi1,deta2,dxi2,dnormeta,dnormxi
    real(dp), dimension(3,3) :: Da,Db
    logical :: bsuccess

    integer :: IVE

    if (iand(celement,int(2**17,I32)) .eq. 0) then

      ! QP1-element with linear mapping to a reference element.

      ! Clear the output array
      !Dbas = 0.0_DP

      ! Calculate the edge midpoints of edges:
      !  DXM(:) := X-coordinates of midpoints
      !  DYM(:) := Y-coordinates of midpoints
      do IVE=1,NVE
        DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,mod(IVE,4)+1))
        DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,mod(IVE,4)+1))
      end do

      ! Calculate the scaling factors for the local coordinate system.
      !  dnormeta := 1 / ||vec_1||_2
      !  dnormxi := 1 / ||vec_2||_2
      !dnormeta = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
      !dnormxi = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)

      ! Calculate the vector eta = (deta1,deta2)
      deta1 = (DXM(2)-DXM(4)) ! * dnormeta
      deta2 = (DYM(2)-DYM(4)) ! * dnormeta

      ! Calculate the vector xi = (dxi1,dxi2)
      dxi1 = (DXM(3)-DXM(1)) ! * dnormxi
      dxi2 = (DYM(3)-DYM(1)) ! * dnormxi

      ! Our basis functions are as follows:
      !
      ! P1(z1,z2) = a1 + b1*z1 + b2*z2 = 1
      ! P2(z1,z2) = a1 + b1*z1 + b2*z2 = z1
      ! P3(z1,z2) = a1 + b1*z1 + b2*z2 = z2
      !
      ! with (z1,z2) the transformed (x,y) in the new coordinate system.
      ! The Pi are defined on the reference element and a linear
      ! mapping sigma:[0,1]^2->R^2 is used to map all the midpoints
      ! from the reference element to the real element:
      !
      !  sigma(0,0) = m
      !  sigma(1,0) = m2 = m + eta/2
      !  sigma(0,1) = m3 = m + xi/2
      !  sigma(-1,0) = m4 = m - eta/2
      !  sigma(0,-1) = m1 = m - xi/2
      !
      !         ^                                       ^
      !   +-----X-----+                        +--------m3--------+
      !   |     |     |                       /         |          \
      !   |     |     |       sigma      ___ /          |           \
      ! --X-----+-----X-->   ------->       m4--------__|m____       \
      !   |     |     |                    /             |    --------m2->
      !   |     |     |                   /              |             \
      !   +-----X-----+                  /               |              \
      !         |                       +----_____       |               \
      !                                           -----__m1_              \
      !                                                  |   -----_____    \
      !                                                                -----+
      !
      ! The basis functions on the real element are defined as
      !
      !  Phi_i(x,y) = Pi(sigma^-1(x,y))
      !
      ! Because sigma is linear, we can calculate the inverse by hand.
      ! sigma is given by the formula
      !
      !   sigma(z1,z2) = m + 1/2 (eta1 xi1) (z1) = m + 1/2 (deta1 dxi1) (z1)
      !                          (eta2 xi2) (z2)           (deta2 dxi2) (z2)
      !
      ! so the inverse mapping is
      !
      !  sigma^-1(z1,z2) = [1/2*(deta1 dxi1)]^-1 * (x - xm)
      !                    [    (deta2 dxi2)]      (y   ym)
      !
      !                  = 2/lambda ( dxi2*(x-xm)  - dxi1*(y-ym)  )
      !                             ( deta1*(y-ym) - deta2*(x-xm) )
      !
      !  with lambda = determinant of the matrix.
      !
      ! So in the first step, calculate (x-xm) and (y-ym).
      dx = (Dpoint(1)-0.5_DP*(DXM(1)+DXM(3)))
      dy = (Dpoint(2)-0.5_DP*(DYM(1)+DYM(3)))

      ! Now calculate the (inverse of the) Jacobian determinant of the linear mapping.
      dxj = 2.0_DP/(deta1*dxi2 - deta2*dxi1)

      ! and calculate the corresponding (eta,xi) = sigma^-1(x,y) of that.
      ! eta is then our first basis function, xi our 2nd:
      !
      !  Phi1(x,y) = 1
      !  Phi2(x,y) = 2/lambda ( dxi2(x-xm) - deta2(y-ym) )
      !  Phi3(x,y) = 2/lambda ( deta1(y-ym) - dxi1(x-xm) )

      Dbas(1,DER_FUNC) = 1
      Dbas(2,DER_FUNC) = dxj * (dxi2*dx - dxi1*dy)
      Dbas(3,DER_FUNC) = dxj * (deta1*dy - deta2*dx)

      if (ubound(Dbas,2) .ge. DER_DERIV_Y) then
        ! Using these definitions, we can also calculate the derivative directly.
        !
        ! X-derivative
        Dbas(1,DER_DERIV_X) = 0.0_DP
        Dbas(2,DER_DERIV_X) = dxj * dxi2
        Dbas(3,DER_DERIV_X) = -dxj * dxi1

        ! Y-derivative
        Dbas(1,DER_DERIV_Y) = 0.0_DP
        Dbas(2,DER_DERIV_Y) = -dxj * deta2
        Dbas(3,DER_DERIV_Y) = dxj * deta1
      end if

    else

      ! QP1-element defining the basis functions directly on the
      ! real element.

      ! Calculate the edge midpoints:
      !  DXM(:) := X-coordinates of midpoints
      !  DYM(:) := Y-coordinates of midpoints
      do IVE=1,NVE
        DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,mod(IVE,4)+1))
        DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,mod(IVE,4)+1))
      end do

      ! Calculate the scaling factors for the local coordinate system.
      !  dnormeta := 1 / ||vec_1||_2
      !  dnormxi := 1 / ||vec_1||_2
      dnormeta = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
      dnormxi = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)

      ! Calculate the vector eta = (deta1,deta2)
      deta1 = (DXM(2)-DXM(4)) * dnormeta
      deta2 = (DYM(2)-DYM(4)) * dnormeta

      ! Calculate the vector xi = (dxi1,dxi2)
      dxi1 = (DXM(3)-DXM(1)) * dnormxi
      dxi2 = (DYM(3)-DYM(1)) * dnormxi

      ! eta/xi define the new coordinate system.
      ! An arbitrary point (x,y) is projected to this coordinate
      ! system using scalar products:
      !
      !  (z1) = sigma(x1,x2) = ( <eta,(x1,x2)> ) = ( eta_1*x1 + eta_2*x2 )
      !  (z2)                  ( <xi, (x1,x2)> )   ( xi_1*x1  + xi_2*x2  )
      !
      !                      = ( eta_1 eta_2 ) ( x1 )
      !                        ( xi_1  xi_2  ) ( x2 )
      !
      ! A general linear polynomial in (eta,xi) has the following form:
      !
      !  p(z1,z2) = a + b*z1 + c*z2
      !
      ! The three local basis functions m1, m2 and m3 of the QP1 approach
      ! are defined by the following restrictions:
      !
      !  m1(0,0) = 1, m1_z1(0,0) = 0, m1_z2(0,0) = 0
      !  m2(0,0) = 0, m2_z1(0,0) = 1, m2_z2(0,0) = 0
      !  m3(0,0) = 0, m3_z1(0,0) = 0, m3_z2(0,0) = 1
      !
      ! (with _z1 and _z2 defining the 1st derivative in direction z1 or
      ! z2,resp.)
      ! Necessarily, our local basis functions in (xi,eta) are defined as follows:
      !
      !  m1(z1,z2) = 1
      !  m2(z1,z2) = z1
      !  m3(z1,z2) = z2
      !
      ! By concatenation with sigma, one receives the basis polynomials
      ! on the real element:
      !
      !  F1(x1,x2) = m1(sigma(x1,x2)) = 1
      !  F2(x1,x2) = m2(sigma(x1,x2)) = eta_1*x1 + eta_2*x2
      !  F3(x1,x2) = m3(sigma(x1,x2)) = xi_1*x1  + xi_2*x2
      !
      ! Each of the three basis functions on the real element is a linear
      ! combination of these three basis polynomials:
      !
      !  Pi(x1,x2) = ai F1(x1,y1)  +  bi F2(x1,x2)  +  ci F3(x1,x2)
      !            = ai  +  bi (eta_1*x1 + eta_2*x2)  +  ci (xi_1*x1  + xi_2*x2)
      !
      ! i.e.
      !
      !  P1(x1,x2) = a1 F1(x1,y1)  +  b1 F2(x1,x2)  +  c1 F3(x1,x2)
      !  P2(x1,x2) = a2 F1(x1,y1)  +  b2 F2(x1,x2)  +  c2 F3(x1,x2)
      !  P3(x1,x2) = a3 F1(x1,y1)  +  b3 F2(x1,x2)  +  c3 F3(x1,x2)
      !
      ! Here Pi satisfies the following restrictions in the midpoint
      ! (xm,ym) of the element:
      !
      !  P1(xm) = 1, P1|eta(xm) = 0, P1|xi(xm) = 0
      !  P2(xm) = 0, P2|eta(xm) = 1, P2|xi(xm) = 0
      !  P3(xm) = 0, P3|eta(xm) = 0, P3|xi(xm) = 1
      !
      ! We need the directional derivative of Pi which is luckily constant:
      !
      !  DPi(x,y) = ( ai F1_x1(x1,y1)  +  bi F2_x1(x1,x2)  +  ci F3_x1(x1,x2) )
      !             ( ai F1_x2(x1,y1)  +  bi F2_x2(x1,x2)  +  ci F3_x2(x1,x2) )
      !
      !           = ( bi eta_1 + ci xi_1 )
      !             ( bi eta_2 + ci xi_2 )
      !
      !  DPi * eta = bi eta_1^2  +  ci xi_1 eta_1  +  bi eta_2^2     +  ci xi_2 eta_2
      !            = bi (eta_1^2 + eta_2^2)        +  ci (eta_1 xi_1 + eta_2 xi_2 )
      !
      !  DPi * xi  = bi eta_1 xi_1  +  ci xi_1^2   +  bi eta_2 xi_2  +  ci xi_2^2
      !            = bi (eta_1 xi_1 + eta_2 xi_2)  +  ci (xi_1^2 + xi_2^2)
      !
      ! The restrictions give us a linear system for the ai, bi and ci:
      !
      !  ( 1     (eta_1*x1 + eta_2*x2)      (xi_1*x1  + xi_2*x2) ) ( a1 a2 a3 ) = ( 1  0  0 )
      !  ( 0       (eta_1^2 + eta_2^2) (eta_1 xi_1 + eta_2 xi_2) ) ( b1 b2 b3 )   ( 0  1  0 )
      !  ( 0 (eta_1 xi_1 + eta_2 xi_2)         (xi_1^2 + xi_2^2) ) ( c1 c2 c3 )   ( 0  0  1 )
      !
      !   ^^^^^^^^^^^^^^^^^^^^^^^ =: A ^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^ =:B ^^^
      !
      ! Calculate the representation of the element midpoint by scalar products
      ! and set up the linear system:

      dmx = 0.5_DP*(DXM(4)+DXM(2))
      dmy = 0.5_DP*(DXM(1)+DXM(3))
      dx = deta1*dmx + deta2*dmy
      dy = dxi1*dmx + dxi2*dmy

      Da(1,1) = 1.0_DP
      Da(2,1) = 0.0_DP
      Da(3,1) = 0.0_DP

      Da(1,2) = deta1*dx + deta2*dy
      Da(2,2) = deta1**2 + deta2**2
      Da(3,2) = deta1*dxi1 + deta2*dxi2

      Da(1,3) = dxi1*dx + dxi2*dy
      Da(2,3) = Da(3,2)
      Da(3,3) = dxi1**2 + dxi2**2

      ! Solve the system to get the coefficients.
      call mprim_invert3x3MatrixDirect(Da,Db,bsuccess)

      ! Use the calculated coefficients to calculate the values of the
      ! basis functions. Take the the point and transform it into
      ! the coordinate system specified by eta and xi:

      dx = deta1*Dpoint(1) + deta2*Dpoint(2)
      dy = dxi1*Dpoint(1) + dxi2*Dpoint(2)

      ! Ok, 1st basis function is clear -- it is constant (as eta/xi are
      ! assumed to be linearly independent and the directional derivatives
      ! are zero). The others have to be calculated using the above coefficients
      ! in Db.
      Dbas(1,DER_FUNC) = 1.0_DP
      Dbas(2,DER_FUNC) = Db(2,2)*dx + Db(3,2)*dy ! + Db(1,2), but this is zero
      Dbas(3,DER_FUNC) = Db(2,3)*dx + Db(3,3)*dy ! + Db(1,3), but this is zero

      if (ubound(Dbas,2) .ge. DER_DERIV_Y) then
        ! X-derivative
        Dbas(1,DER_DERIV_X) = 0.0_DP
        Dbas(2,DER_DERIV_X) = Db(2,2)*deta1 + Db(3,2)*dxi1
        Dbas(3,DER_DERIV_X) = Db(2,3)*deta1 + Db(3,3)*dxi1

        ! Y-derivative
        Dbas(1,DER_DERIV_Y) = 0.0_DP
        Dbas(2,DER_DERIV_Y) = Db(2,2)*deta2 + Db(3,2)*dxi2
        Dbas(3,DER_DERIV_Y) = Db(2,3)*deta2 + Db(3,3)*dxi2
      end if

    end if

  end subroutine

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
!   |           |                       /                    \           .
!   |           |       sigma          /                      \          .
!   e4          e2     ------->       E4                       \         .
!   |           |                    /                          E2
!   |           |                   /                            \       .
!   +-----e1----+                  /                              \      .
!                                 +----_____                       \     .
!                                           -----__E1_              \    .
!                                                      -----_____    \   .
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
!          /          ^         \                                       .
!         /           |vec_2     \                                      .
!        X--------____|___        \                                     .
!       /             |   -------->X                                    .
!      /              |     vec_1   \                                   .
!     /               |              \                                  .
!    +----_____       |               \                                 .
!              -----__X__              \                                .
!                         -----_____    \                               .
!                                   -----+
!
! We shift both of these vectors vec_1 and vec_2 into the origin (0,0)
! and normalise them to have length=1 (otherwise, the LBB condition might
! get violated!). So the picture we are looking at here is the following:
!
!   ^ xi             +---------X--------+
!   |               /          ^         \                               .
!   |              /           |vec_2     \                              .
!   |             X--------____|___        \                             .
!   |            /             |   -------->X                            .
!   |           /              |     vec_1   \                           .
!   |          /               |              \                          .
!   |         +----_____       |               \                         .
!   |                   -----__X__              \                        .
!   |                              -----_____    \                       .
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
! Now remember, that eta/xi has length = 1. Then, we can get the new
! coordinates of (x,y) by using scalar products:
!
!   z1 = < eta, (x,y) >
!   z2 = < xi, (x,y) >
!
! This is a linear (not bilinear!) mapping of (x,y) into
! the new coordinate system.
!
! So we end up with:
!
!   ( z1 ) := r(x,y) := ( eta_1 eta_2 ) ( x )
!   ( z2 )              ( xi_1  xi_2  ) ( y )
!
!                     =  ( eta_1 x  +  eta_2 y )
!                        ( xi_1 x   +  xi_2 y )
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
! To evaluate these mi`s in the new coordinate system, we concatenate them
! with the mapping r(.,.). As result, we get the four functions F1,F2,F3,F4
! which are defined as functions at the bottom of this routine:
!
!  F1(x,y) := m1(r(x,y)) = 1
!  F2(x,y) := m2(r(x,y)) = eta_1 x + eta_2 y
!  F3(x,y) := m3(r(x,y)) = xi_1 x  + xi_2 y
!  F4(x,y) := m4(r(x,y)) = ( eta_1 x  +  eta_2 y )^2 - ( xi_1 x  +  xi_2 y )^2
!                        =            ( eta_1^2 - xi_1^2 ) x^2
!                          + 2 ( eta_1 eta_2 - xi_1 xi_2 ) x y
!                          +          ( eta_2^2 - xi_2^2 ) y^2
!
! So the polynomials have now the form:
!
!  P1(r(x,y)) = a1 F1(x,y)  +  b1 F2(x,y)  +  c1 F3(x,y)  +  d1 F4(x,y)
!  P2(r(x,y)) = a2 F1(x,y)  +  b2 F2(x,y)  +  c2 F3(x,y)  +  d2 F4(x,y)
!  P3(r(x,y)) = a3 F1(x,y)  +  b3 F2(x,y)  +  c3 F3(x,y)  +  d3 F4(x,y)
!  P4(r(x,y)) = a4 F1(x,y)  +  b4 F2(x,y)  +  c4 F3(x,y)  +  d4 F4(x,y)
!
! It does not matter whether the local coordinate system starts in (0,0) or in
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
!               ^^^^^ =:B ^^^^^
!
! So to get all the coefficients, one has to calculate B = V^-1 !
! The entries of the matrix V = {v_ij} are defined (because of the linearity
! of the integral) as
!
!         vij = 1/|ei| int_ei Fj(x,y) ds
!
! Now let us go...

!**************************************************************************
! Element subroutines for nonparametric Q1~ element, integral mean value
! based.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>

  pure subroutine elem_EM30 (celement, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE),
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(in) :: ddetj


  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on the real element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables
  integer :: IVE,IA,IK
  real(DP) :: PXL,PYL,PXU,PYU, D1,D2
  real(DP),dimension(4) :: DXM,DYM,DLX,DLY
  real(DP),dimension(4,4) :: A       ! local 4x4 system
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  real(DP),dimension(4,6) :: COB
  real(DP) :: dx, dy
  logical :: bsuccess

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

  do IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,mod(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,mod(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1)-Dcoords(2,IVE))
  end do

  ! Calculate the scaling factors for the local coordinate system.
  !  D1 := 1 / ||xi||_2
  !  D2 := 1 / ||eta||_2

  D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)

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
  do IA = 1,NVE

    ! Calculate the X- and Y-coordinates of the two Gauss points on the
    ! current edge. (PXL,PYL) = 1st Gauss point, (PXU,PYU) = 2nd Gauss point.

    PXL = DXM(IA)-sqrt(1.0_DP/3.0_DP)*DLX(IA)
    PYL = DYM(IA)-sqrt(1.0_DP/3.0_DP)*DLY(IA)
    PXU = DXM(IA)+sqrt(1.0_DP/3.0_DP)*DLX(IA)
    PYU = DYM(IA)+sqrt(1.0_DP/3.0_DP)*DLY(IA)

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
  end do

  ! Invert that matrix V to get the matrix of the coefficients of the
  ! four polynomials. The matix A (=V) is replaced by the inverse.
  call mprim_invertMatrixPivot(A,4,bsuccess)

  ! Ok, the coefficients ai, bi, ci, di of the matrix B are calculated.
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

  do IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
  end do

  ! Function values

  if (BDER(DER_FUNC)) then
    Dbas(1,DER_FUNC)= COB(1,1)*dx**2+COB(1,2)*dy**2+COB(1,3)*dx*dy &
                     +COB(1,4)*dx   +COB(1,5)*dy   +COB(1,6)
    Dbas(2,DER_FUNC)= COB(2,1)*dx**2+COB(2,2)*dy**2+COB(2,3)*dx*dy &
                     +COB(2,4)*dx   +COB(2,5)*dy   +COB(2,6)
    Dbas(3,DER_FUNC)= COB(3,1)*dx**2+COB(3,2)*dy**2+COB(3,3)*dx*dy &
                     +COB(3,4)*dx   +COB(3,5)*dy   +COB(3,6)
    Dbas(4,DER_FUNC)= COB(4,1)*dx**2+COB(4,2)*dy**2+COB(4,3)*dx*dy &
                     +COB(4,4)*dx   +COB(4,5)*dy   +COB(4,6)
  end if

  ! Derivatives:

  if (BDER(DER_DERIV_X) .or. BDER(DER_DERIV_Y)) then
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
  end if

  contains

    ! Auxiliary functions

    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function

    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function

    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function

    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_EM30_mult (celement, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given points. The coordinates are expected
  ! on the real element!
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
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
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables
  logical :: bsuccess
  integer :: IVE,IA,IK,i
  real(DP) :: PXL,PYL,PXU,PYU, D1,D2
  real(DP),dimension(4) :: DXM,DYM,DLX,DLY
  real(DP),dimension(4,4) :: A       ! local 4x4 system
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  real(DP),dimension(4,6) :: COB

  ! Clear the output array
  !Dbas = 0.0_DP

  do IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,mod(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,mod(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1)-Dcoords(2,IVE))
  end do

  D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
  CA1 = (DXM(2)-DXM(4)) * D1
  CB1 = (DYM(2)-DYM(4)) * D1
  CA2 = (DXM(3)-DXM(1)) * D2
  CB2 = (DYM(3)-DYM(1)) * D2
  CA3 = CA1**2-CA2**2
  CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
  CC3 = CB1**2-CB2**2

  do IA = 1,4
    PXL = DXM(IA)-sqrt(1.0_DP/3.0_DP)*DLX(IA)
    PYL = DYM(IA)-sqrt(1.0_DP/3.0_DP)*DLY(IA)
    PXU = DXM(IA)+sqrt(1.0_DP/3.0_DP)*DLX(IA)
    PYU = DYM(IA)+sqrt(1.0_DP/3.0_DP)*DLY(IA)
    A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
  end do

  call mprim_invertMatrixPivot(A,4,bsuccess)

  do IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
  end do

  ! Function values

  if (BDER(DER_FUNC)) then
    do i=1,npoints
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
    end do
  end if

  ! Derivatives:

  if (BDER(DER_DERIV_X) .or. BDER(DER_DERIV_Y)) then
    ! x-derivatives

    do i=1,npoints
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
    end do
  end if

  contains

    ! Auxiliary functions

    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function

    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function

    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function

    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine elem_EM30_sim (celement, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, &
                            Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
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
  ! The coordinates are expected on the real element.
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th
  !   basis function of the finite element in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables
  logical :: bsuccess
  integer :: IVE,IA,IK,i,j
  real(DP) :: PXL,PYL,PXU,PYU, D1,D2
  real(DP),dimension(4) :: DXM,DYM,DLX,DLY
  real(DP), dimension(4,4) :: A       ! local 4x4 system
  real(DP), dimension(4,4) :: CK
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  real(DP) :: COB(6,4)
  integer :: npointsfunc,npointsderx,npointsdery

  ! Calculate the loop counters in advance. Help us to get rid
  ! of any if-commands in the element loop below.
  npointsfunc = 0
  npointsderx = 0
  npointsdery = 0
  if (Bder(DER_FUNC)) npointsfunc = npoints
  if (Bder(DER_DERIV_X)) npointsderx = npoints
  if (Bder(DER_DERIV_Y)) npointsdery = npoints

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Check which element variant we have; with or without pivoting...
  if (iand(celement,int(2**17,I32)) .eq. 0) then

    ! Check whether to scale the local coordinate system or not.

    if (iand(celement,int(2**18,I32)) .eq. 0) then

      ! Use pivoting and scaled local coordinate system for
      ! increased numerical stability.

      ! Loop over the elements
      !$omp parallel do default(shared)&
      !$omp private(A,CA1,CA2,CA3,CB1,CB2,CB3,CC3,COB,D1,D2,DLX,DLY,DXM,DYM,&
      !$omp         IA,IK,IVE,PXL,PXU,PYL,PYU,bsuccess,i)&
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
        end do

        D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
        D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
        CA1 = (DXM(2)-DXM(4)) * D1
        CB1 = (DYM(2)-DYM(4)) * D1
        CA2 = (DXM(3)-DXM(1)) * D2
        CB2 = (DYM(3)-DYM(1)) * D2
        CA3 = CA1**2-CA2**2
        CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
        CC3 = CB1**2-CB2**2

        do IA = 1,4
          PXL = DXM(IA)-sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYL = DYM(IA)-sqrt(1.0_DP/3.0_DP)*DLY(IA)
          PXU = DXM(IA)+sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYU = DYM(IA)+sqrt(1.0_DP/3.0_DP)*DLY(IA)

          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        end do

        ! Invert the matrix in-place.
        call mprim_invertMatrixPivot(A,4,bsuccess)

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

        do IK=1,4
          COB(1,IK) = A(IK,4)*CA3
          COB(2,IK) = A(IK,4)*CC3
          COB(3,IK) = A(IK,4)*CB3
          COB(4,IK) = A(IK,2)*CA1+A(IK,3)*CA2
          COB(5,IK) = A(IK,2)*CB1+A(IK,3)*CB2
          COB(6,IK) = A(IK,1)
        end do

        ! Function values
        do i=1,npointsfunc   ! either 0 or npoints

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
        end do

        ! x-derivatives

        do i=1,npointsderx   ! either 0 or npoints
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
        end do

        ! y-derivatives

        do i=1,npointsdery   ! either 0 or npoints
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

        end do

      end do ! ielement
      !$omp end parallel do

    else

      ! Use pivoting for increased numerical stability.
      ! Do not scaled local coordinate system

      ! Loop over the elements
      !$omp parallel do default(shared)&
      !$omp private(A,CA1,CA2,CA3,CB1,CB2,CB3,CC3,COB,DLX,DLY,DXM,DYM,&
      !$omp         IA,IK,IVE,PXL,PXU,PYL,PYU,bsuccess,i)&
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
        end do

        ! Do not scale the local coordinate system in this approach.
        ! Saves a little it numerical effort.

        CA1 = (DXM(2)-DXM(4))
        CB1 = (DYM(2)-DYM(4))
        CA2 = (DXM(3)-DXM(1))
        CB2 = (DYM(3)-DYM(1))
        CA3 = CA1**2-CA2**2
        CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
        CC3 = CB1**2-CB2**2

        do IA = 1,4
          PXL = DXM(IA)-sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYL = DYM(IA)-sqrt(1.0_DP/3.0_DP)*DLY(IA)
          PXU = DXM(IA)+sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYU = DYM(IA)+sqrt(1.0_DP/3.0_DP)*DLY(IA)

          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        end do

        ! Invert the matrix in-place.
        call mprim_invertMatrixPivot(A,4,bsuccess)

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

        do IK=1,4
          COB(1,IK) = A(IK,4)*CA3
          COB(2,IK) = A(IK,4)*CC3
          COB(3,IK) = A(IK,4)*CB3
          COB(4,IK) = A(IK,2)*CA1+A(IK,3)*CA2
          COB(5,IK) = A(IK,2)*CB1+A(IK,3)*CB2
          COB(6,IK) = A(IK,1)
        end do

        ! Function values
        do i=1,npointsfunc   ! either 0 or npoints

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
        end do

        ! x-derivatives

        do i=1,npointsderx   ! either 0 or npoints
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
        end do

        ! y-derivatives

        do i=1,npointsdery   ! either 0 or npoints
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

        end do

      end do ! ielement
      !$omp end parallel do

    end if

  else

    ! Do not use pivoting.

    ! Check whether to scae the local coordinate system or not.

    if (iand(celement,int(2**18,I32)) .eq. 0) then

      ! Loop over the elements
      !$omp parallel do default(shared)&
      !$omp private(A,CA1,CA2,CA3,CB1,CB2,CB3,CC3,CK,COB,D1,D2,DLX,DLY,DXM,DYM,&
      !$omp         IA,IK,IVE,PXL,PXU,PYL,PYU,bsuccess,i)&
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
        end do

        D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
        D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
        CA1 = (DXM(2)-DXM(4)) * D1
        CB1 = (DYM(2)-DYM(4)) * D1
        CA2 = (DXM(3)-DXM(1)) * D2
        CB2 = (DYM(3)-DYM(1)) * D2
        CA3 = CA1**2-CA2**2
        CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
        CC3 = CB1**2-CB2**2

        do IA = 1,4
          PXL = DXM(IA)-sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYL = DYM(IA)-sqrt(1.0_DP/3.0_DP)*DLY(IA)
          PXU = DXM(IA)+sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYU = DYM(IA)+sqrt(1.0_DP/3.0_DP)*DLY(IA)
          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        end do

        ! Invert the matrix to get the coefficients.
        ! Use direct inversion and save the result to CK directly.
        call mprim_invert4x4MatrixDirect(A,CK,bsuccess)

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
        do IK=1,4
          COB(1,IK) = CK(IK,4)*CA3
          COB(2,IK) = CK(IK,4)*CC3
          COB(3,IK) = CK(IK,4)*CB3
          COB(4,IK) = CK(IK,2)*CA1+CK(IK,3)*CA2
          COB(5,IK) = CK(IK,2)*CB1+CK(IK,3)*CB2
          COB(6,IK) = CK(IK,1)
        end do

        ! Function values
        do i=1,npointsfunc   ! either 0 or npoints

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
        end do

        ! x-derivatives

        do i=1,npointsderx   ! either 0 or npoints
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
        end do

        ! y-derivatives

        do i=1,npointsdery   ! either 0 or npoints
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

        end do

      end do ! ielement
      !$omp end parallel do

    else

      ! Loop over the elements
      !$omp parallel do default(shared)&
      !$omp private(A,CA1,CA2,CA3,CB1,CB2,CB3,CC3,CK,COB,DLX,DLY,DXM,DYM,&
      !$omp         IA,IVE,PXL,PXU,PYL,PYU,bsuccess,i)&
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
        end do

        ! Do not scale the local coordinate in this approach.
        ! Saves a little bit numerical effort.

        CA1 = (DXM(2)-DXM(4))
        CB1 = (DYM(2)-DYM(4))
        CA2 = (DXM(3)-DXM(1))
        CB2 = (DYM(3)-DYM(1))
        CA3 = CA1**2-CA2**2
        CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
        CC3 = CB1**2-CB2**2

        do IA = 1,4
          PXL = DXM(IA)-sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYL = DYM(IA)-sqrt(1.0_DP/3.0_DP)*DLY(IA)
          PXU = DXM(IA)+sqrt(1.0_DP/3.0_DP)*DLX(IA)
          PYU = DYM(IA)+sqrt(1.0_DP/3.0_DP)*DLY(IA)
          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        end do

        ! Invert the matrix to get the coefficients.
        ! Use direct inversion and save the result to CK directly.
        call mprim_invert4x4MatrixDirect(A,CK,bsuccess)

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
        do IK=1,4
          COB(1,IK) = CK(IK,4)*CA3
          COB(2,IK) = CK(IK,4)*CC3
          COB(3,IK) = CK(IK,4)*CB3
          COB(4,IK) = CK(IK,2)*CA1+CK(IK,3)*CA2
          COB(5,IK) = CK(IK,2)*CB1+CK(IK,3)*CB2
          COB(6,IK) = CK(IK,1)
        end do

        ! Function values
        do i=1,npointsfunc   ! either 0 or npoints

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
        end do

        ! x-derivatives

        do i=1,npointsderx   ! either 0 or npoints
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
        end do

        ! y-derivatives

        do i=1,npointsdery   ! either 0 or npoints
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

        end do

      end do ! ielement
      !$omp end parallel do

    end if

  end if

  contains

    ! Auxiliary functions

    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function

    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function

    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function

    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine

!**************************************************************************
! Element subroutines for parametric Q1~ element, integral mean value
! based.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>

  pure subroutine elem_E030 (celement, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D) :: Dhelp
  real(DP) :: dx,dy

  real(DP) :: dxj !auxiliary variable

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
  ! That is even faster than when using three IF commands for preventing
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

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_E030_mult (celement, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC,i) = 0.125E0_DP* &
          (-3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)-4.0_DP*Dpoints(2,i)+2.0_DP)
      Dbas(2,DER_FUNC,i) = 0.125E0_DP* &
          ( 3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)+4.0_DP*Dpoints(1,i)+2.0_DP)
      Dbas(3,DER_FUNC,i) = -0.125E0_DP* &
          ( 3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)-4.0_DP*Dpoints(2,i)-2.0_DP)
      Dbas(4,DER_FUNC,i) = -0.125E0_DP* &
          (-3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)+4.0_DP*Dpoints(1,i)-2.0_DP)
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj(:) = 0.125E0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
      Dhelp(1,1,i) = -6.0_DP*Dpoints(1,i)
      Dhelp(2,1,i) =  6.0_DP*Dpoints(1,i)+4.0_DP
      Dhelp(3,1,i) = -6.0_DP*Dpoints(1,i)
      Dhelp(4,1,i) =  6.0_DP*Dpoints(1,i)-4.0_DP
      Dhelp(1,2,i) =  6.0_DP*Dpoints(2,i)-4.0_DP
      Dhelp(2,2,i) = -6.0_DP*Dpoints(2,i)
      Dhelp(3,2,i) =  6.0_DP*Dpoints(2,i)+4.0_DP
      Dhelp(4,2,i) = -6.0_DP*Dpoints(2,i)
    end do

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_E030_sim (celement, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, &
                            Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 0.125E0_DP* &
            (-3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-4.0_DP*Dpoints(2,i,j)+2.0_DP)
        Dbas(2,DER_FUNC,i,j) = 0.125E0_DP* &
            ( 3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+4.0_DP*Dpoints(1,i,j)+2.0_DP)
        Dbas(3,DER_FUNC,i,j) = -0.125E0_DP* &
            ( 3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-4.0_DP*Dpoints(2,i,j)-2.0_DP)
        Dbas(4,DER_FUNC,i,j) = -0.125E0_DP* &
            (-3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+4.0_DP*Dpoints(1,i,j)-2.0_DP)
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj,Dhelp) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 0.125E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
        Dhelp(1,1,i) = -6.0_DP*Dpoints(1,i,j)
        Dhelp(2,1,i) =  6.0_DP*Dpoints(1,i,j)+4.0_DP
        Dhelp(3,1,i) = -6.0_DP*Dpoints(1,i,j)
        Dhelp(4,1,i) =  6.0_DP*Dpoints(1,i,j)-4.0_DP
        Dhelp(1,2,i) =  6.0_DP*Dpoints(2,i,j)-4.0_DP
        Dhelp(2,2,i) = -6.0_DP*Dpoints(2,i,j)
        Dhelp(3,2,i) =  6.0_DP*Dpoints(2,i,j)+4.0_DP
        Dhelp(4,2,i) = -6.0_DP*Dpoints(2,i,j)
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(1,1,i) &
                                    - Djac(2,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(2,1,i) &
                                    - Djac(2,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(3,1,i) &
                                    - Djac(2,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(4,1,i) &
                                    - Djac(2,i,j) * Dhelp(4,2,i))
        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(1,1,i) &
                                    - Djac(1,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(2,1,i) &
                                    - Djac(1,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(3,1,i) &
                                    - Djac(1,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(4,1,i) &
                                    - Djac(1,i,j) * Dhelp(4,2,i))
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric Q1~ element with bubble, integral mean
! value based.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>

  pure subroutine elem_EB30 (celement, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1TB.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(5,NDIM2D) :: Dhelp
  real(DP) :: dx,dy,dxy

  real(DP) :: dxj !auxiliary variable

  ! The Q1 element with bubble is specified by 5 polynomials on
  ! the reference element.
  ! These 5 polynomials are:
  !
  !  p_1(x,y) = -3/8 (x^2-y^2)  -  1/2 y  +  1/4
  !  p_2(x,y) =  3/8 (x^2-y^2)  +  1/2 x  +  1/4
  !  p_3(x,y) = -3/8 (x^2-y^2)  +  1/2 y  +  1/4
  !  p_4(x,y) =  3/8 (x^2-y^2)  -  1/2 x  +  1/4
  !  p_5(x,y) =  9 xy
  !
  ! These five polynomials are constructed in a way such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,5
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (p_i(Ej(t))) d(t) = kronecker(i,j  ) * |ej|
  !   }
  !   Int_T (p_i(x,y)*x*y) d(x,y) = kronecker(i,5) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th edge of the reference quadrilateral
  ! Ej: [-1,1] -> ej being the parametrisation of an edge ej
  ! |ej| = 2 being the length of the edge ej
  ! T being the reference quadrilateral
  ! |T| = 4 being the area of the reference quadrilateral

  ! Clear the output array
  !Dbas = 0.0_DP

  dx = Dpoint(1)
  dy = Dpoint(2)
  dxy = dx**2 - dy**2

  ! Remark: The Q1~-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) =  0.125_DP*(-3.0_DP*dxy - 4.0_DP*dy + 2.0_DP)
    Dbas(2,DER_FUNC) =  0.125_DP*( 3.0_DP*dxy + 4.0_DP*dx + 2.0_DP)
    Dbas(3,DER_FUNC) = -0.125_DP*( 3.0_DP*dxy - 4.0_DP*dy - 2.0_DP)
    Dbas(4,DER_FUNC) = -0.125_DP*(-3.0_DP*dxy + 4.0_DP*dx - 2.0_DP)
    Dbas(5,DER_FUNC) =  9.0_DP * dx * dy
!  endif

  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 0.125E0_DP / ddetj

    ! Remember that we defined dxj = (1/8) / det instead of 1/det, therefore
    ! we have to multiply the derivatives of the fifth basis function by 8 !

    ! x- and y-derivatives on reference element
    Dhelp(1,1) = -6.0_DP*dx
    Dhelp(2,1) =  6.0_DP*dx + 4.0_DP
    Dhelp(3,1) = -6.0_DP*dx
    Dhelp(4,1) =  6.0_DP*dx - 4.0_DP
    Dhelp(5,1) = 72.0_DP*dy

    Dhelp(1,2) =  6.0_DP*dy - 4.0_DP
    Dhelp(2,2) = -6.0_DP*dy
    Dhelp(3,2) =  6.0_DP*dy + 4.0_DP
    Dhelp(4,2) = -6.0_DP*dy
    Dhelp(5,2) = 72.0_DP*dx

    ! x-derivatives on current element
!    if (Bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(1,1) - Djac(2) * Dhelp(1,2))
      Dbas(2,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(2,1) - Djac(2) * Dhelp(2,2))
      Dbas(3,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(3,1) - Djac(2) * Dhelp(3,2))
      Dbas(4,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(4,1) - Djac(2) * Dhelp(4,2))
      Dbas(5,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(5,1) - Djac(2) * Dhelp(5,2))
!    endif

    ! y-derivatives on current element
!    if (Bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(1,1) - Djac(1) * Dhelp(1,2))
      Dbas(2,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(2,1) - Djac(1) * Dhelp(2,2))
      Dbas(3,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(3,1) - Djac(1) * Dhelp(3,2))
      Dbas(4,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(4,1) - Djac(1) * Dhelp(4,2))
      Dbas(5,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(5,1) - Djac(1) * Dhelp(5,2))
!    endif
!  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_EB30_mult (celement, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1TB.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(5,NDIM2D,npoints) :: Dhelp
  real(DP) :: dx,dy,dxy

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints

      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      dxy = dx**2 - dy**2

      Dbas(1,DER_FUNC,i) =  0.125_DP*(-3.0_DP*dxy - 4.0_DP*dy + 2.0_DP)
      Dbas(2,DER_FUNC,i) =  0.125_DP*( 3.0_DP*dxy + 4.0_DP*dx + 2.0_DP)
      Dbas(3,DER_FUNC,i) = -0.125_DP*( 3.0_DP*dxy - 4.0_DP*dy - 2.0_DP)
      Dbas(4,DER_FUNC,i) = -0.125_DP*(-3.0_DP*dxy + 4.0_DP*dx - 2.0_DP)
      Dbas(5,DER_FUNC,i) =  9.0_DP * dx * dy
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj(:) = 0.125E0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dhelp(1,1,i) = -6.0_DP*dx
      Dhelp(2,1,i) =  6.0_DP*dx + 4.0_DP
      Dhelp(3,1,i) = -6.0_DP*dx
      Dhelp(4,1,i) =  6.0_DP*dx - 4.0_DP
      Dhelp(5,1,i) = 72.0_DP*dy
      Dhelp(1,2,i) =  6.0_DP*dy - 4.0_DP
      Dhelp(2,2,i) = -6.0_DP*dy
      Dhelp(3,2,i) =  6.0_DP*dy + 4.0_DP
      Dhelp(4,2,i) = -6.0_DP*dy
      Dhelp(5,2,i) = 72.0_DP*dx
    end do

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        Dbas(1,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(1,1,i) - Djac(2,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(2,1,i) - Djac(2,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(3,1,i) - Djac(2,i) * Dhelp(3,2,i))
        Dbas(4,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(4,1,i) - Djac(2,i) * Dhelp(4,2,i))
        Dbas(5,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(5,1,i) - Djac(2,i) * Dhelp(5,2,i))
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
        Dbas(5,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(5,1,i) - Djac(1,i) * Dhelp(5,2,i))
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_EB30_sim (celement, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, &
                            Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1TB.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(5,NDIM2D,npoints) :: Dhelp
  real(DP) :: dx,dy,dxy

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i,dx,dy,dxy) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dxy = dx**2 - dy**2
        Dbas(1,DER_FUNC,i,j) =  0.125_DP*(-3.0_DP*dxy - 4.0_DP*dy + 2.0_DP)
        Dbas(2,DER_FUNC,i,j) =  0.125_DP*( 3.0_DP*dxy + 4.0_DP*dx + 2.0_DP)
        Dbas(3,DER_FUNC,i,j) = -0.125_DP*( 3.0_DP*dxy - 4.0_DP*dy - 2.0_DP)
        Dbas(4,DER_FUNC,i,j) = -0.125_DP*(-3.0_DP*dxy + 4.0_DP*dx - 2.0_DP)
        Dbas(5,DER_FUNC,i,j) =  9.0_DP * dx * dy
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj,dx,dy,Dhelp) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 0.125E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dhelp(1,1,i) = -6.0_DP*dx
        Dhelp(2,1,i) =  6.0_DP*dx + 4.0_DP
        Dhelp(3,1,i) = -6.0_DP*dx
        Dhelp(4,1,i) =  6.0_DP*dx - 4.0_DP
        Dhelp(5,1,i) = 72.0_DP*dy
        Dhelp(1,2,i) =  6.0_DP*dy - 4.0_DP
        Dhelp(2,2,i) = -6.0_DP*dy
        Dhelp(3,2,i) =  6.0_DP*dy + 4.0_DP
        Dhelp(4,2,i) = -6.0_DP*dy
        Dhelp(5,2,i) = 72.0_DP*dx
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(1,1,i) &
                                    - Djac(2,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(2,1,i) &
                                    - Djac(2,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(3,1,i) &
                                    - Djac(2,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(4,1,i) &
                                    - Djac(2,i,j) * Dhelp(4,2,i))
          Dbas(5,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(5,1,i) &
                                    - Djac(2,i,j) * Dhelp(5,2,i))
        !end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        !do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(1,1,i) &
                                    - Djac(1,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(2,1,i) &
                                    - Djac(1,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(3,1,i) &
                                    - Djac(1,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(4,1,i) &
                                    - Djac(1,i,j) * Dhelp(4,2,i))
          Dbas(5,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(5,1,i) &
                                    - Djac(1,i,j) * Dhelp(5,2,i))
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

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

  pure subroutine elem_EM31 (celement, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on the real element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables
  logical :: bsuccess
  integer :: IVE,IA,IK
  real(DP) :: D1,D2
  real(DP),dimension(4) :: DXM,DYM,DLX,DLY
  real(DP),dimension(4,4) :: A       ! local 4x4 system
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  real(DP),dimension(4,6) :: COB
  real(DP) :: dx, dy

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

  do IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,mod(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,mod(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1)-Dcoords(2,IVE))
  end do

  ! Calculate the scaling factors for the local coordinate system.
  !  D1 := 1 / ||xi||_2
  !  D2 := 1 / ||eta||_2

  D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)

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

  ! Calculate the matrix V (=A) with vij = Fj (m_i).
  ! Loop over the edges.
  do IA = 1,NVE
    ! Set up the coefficients of the linear system to calculate ai, bi, ci and di.
    ! Use the X- and Y-coordinates of the midpoint of every edge to evaluate
    ! the Fi.
    A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
  end do

  ! Invert that matrix V to get the matrix of the coefficients of the
  ! four polynomials. The matix A (=V) is replaced by the inverse.
  call mprim_invertMatrixPivot(A,4,bsuccess)

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

  do IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
  end do

  ! Function values

  if (BDER(DER_FUNC)) then
    Dbas(1,DER_FUNC)= COB(1,1)*dx**2+COB(1,2)*dy**2+COB(1,3)*dx*dy &
                     +COB(1,4)*dx   +COB(1,5)*dy   +COB(1,6)
    Dbas(2,DER_FUNC)= COB(2,1)*dx**2+COB(2,2)*dy**2+COB(2,3)*dx*dy &
                     +COB(2,4)*dx   +COB(2,5)*dy   +COB(2,6)
    Dbas(3,DER_FUNC)= COB(3,1)*dx**2+COB(3,2)*dy**2+COB(3,3)*dx*dy &
                     +COB(3,4)*dx   +COB(3,5)*dy   +COB(3,6)
    Dbas(4,DER_FUNC)= COB(4,1)*dx**2+COB(4,2)*dy**2+COB(4,3)*dx*dy &
                     +COB(4,4)*dx   +COB(4,5)*dy   +COB(4,6)
  end if

  ! Derivatives:

  if (BDER(DER_DERIV_X) .or. BDER(DER_DERIV_Y)) then
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
  end if

  contains

    ! Auxiliary functions

    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function

    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function

    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function

    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_EM31_mult (celement, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given points. The coordinates are expected
  ! on the real element!
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
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
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables
  logical :: bsuccess
  integer :: IVE,IA,IK,i
  real(DP) :: D1,D2
  real(DP),dimension(4) :: DXM,DYM,DLX,DLY
  real(DP),dimension(4,4) :: A       ! local 4x4 system
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  real(DP),dimension(4,6) :: COB

  ! Clear the output array
  !Dbas = 0.0_DP

  do IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,mod(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,mod(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1)-Dcoords(2,IVE))
  end do

  D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
  CA1 = (DXM(2)-DXM(4)) * D1
  CB1 = (DYM(2)-DYM(4)) * D1
  CA2 = (DXM(3)-DXM(1)) * D2
  CB2 = (DYM(3)-DYM(1)) * D2
  CA3 = CA1**2-CA2**2
  CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
  CC3 = CB1**2-CB2**2

  do IA = 1,4
    A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
  end do

  call mprim_invertMatrixPivot(A,4,bsuccess)

  do IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
  end do

  ! Function values

  if (BDER(DER_FUNC)) then
    do i=1,npoints
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
    end do
  end if

  ! Derivatives:

  if (BDER(DER_DERIV_X) .or. BDER(DER_DERIV_Y)) then
    ! x-derivatives

    do i=1,npoints
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
    end do
  end if

  contains

    ! Auxiliary functions

    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function

    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function

    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function

    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine elem_EM31_sim (celement, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, &
                            Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
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
  ! The coordinates are expected on the real element.
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th
  !   basis function of the finite element in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables
  logical :: bsuccess
  integer :: IVE,IA,IK,i,j
  real(DP) :: D1,D2
  real(DP),dimension(4) :: DXM,DYM,DLX,DLY
  real(DP), dimension(4,4) :: A       ! local 4x4 system
  real(DP), dimension(4,4) :: CK
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  real(DP) :: COB(6,4)
  integer :: npointsfunc,npointsderx,npointsdery

  ! Calculate the loop counters in advance. Help us to get rid
  ! of any if-commands in the element loop below.
  npointsfunc = 0
  npointsderx = 0
  npointsdery = 0
  if (Bder(DER_FUNC)) npointsfunc = npoints
  if (Bder(DER_DERIV_X)) npointsderx = npoints
  if (Bder(DER_DERIV_Y)) npointsdery = npoints

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Check which element variant we have; with or without pivoting...
  if (iand(celement,int(2**17,I32)) .eq. 0) then

    ! Use pivoting for increased numerical stability.

    ! Loop over the elements
    !$omp parallel do default(shared)&
    !$omp private(A,CA1,CA2,CA3,CB1,CB2,CB3,CC3,COB,D1,D2,DLX,DLY,DXM,DYM,&
    !$omp         IA,IK,IVE,bsuccess,i)&
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do IVE=1,NVE
        DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
        DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
        DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
        DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
      end do

      D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
      D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
      CA1 = (DXM(2)-DXM(4)) * D1
      CB1 = (DYM(2)-DYM(4)) * D1
      CA2 = (DXM(3)-DXM(1)) * D2
      CB2 = (DYM(3)-DYM(1)) * D2
      CA3 = CA1**2-CA2**2
      CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
      CC3 = CB1**2-CB2**2

      do IA = 1,4
        A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
      end do

      ! Invert the matrix in-place.
      call mprim_invertMatrixPivot(A,4,bsuccess)

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

      do IK=1,4
        COB(1,IK) = A(IK,4)*CA3
        COB(2,IK) = A(IK,4)*CC3
        COB(3,IK) = A(IK,4)*CB3
        COB(4,IK) = A(IK,2)*CA1+A(IK,3)*CA2
        COB(5,IK) = A(IK,2)*CB1+A(IK,3)*CB2
        COB(6,IK) = A(IK,1)
      end do

      ! Function values
      do i=1,npointsfunc   ! either 0 or npoints

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
      end do

      ! x-derivatives

      do i=1,npointsderx   ! either 0 or npoints
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
      end do

      ! y-derivatives

      do i=1,npointsdery   ! either 0 or npoints
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

      end do

    end do ! ielement
    !$omp end parallel do

  else

    ! Do not use pivoting.

    ! Loop over the elements
    !$omp parallel do default(shared)&
    !$omp private(A,CA1,CA2,CA3,CB1,CB2,CB3,CC3,CK,COB,D1,D2,DLX,DLY,DXM,DYM,&
    !$omp         IA,IK,IVE,bsuccess,i)&
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do IVE=1,NVE
        DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
        DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
        DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
        DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
      end do

      D1 = 1.0_DP / sqrt((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
      D2 = 1.0_DP / sqrt((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
      CA1 = (DXM(2)-DXM(4)) * D1
      CB1 = (DYM(2)-DYM(4)) * D1
      CA2 = (DXM(3)-DXM(1)) * D2
      CB2 = (DYM(3)-DYM(1)) * D2
      CA3 = CA1**2-CA2**2
      CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
      CC3 = CB1**2-CB2**2

      do IA = 1,4
        A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
      end do

      ! Invert the matrix to get the coefficients.
      ! Use direct inversion and save the result to CK directly.
      call mprim_invert4x4MatrixDirect(A,CK,bsuccess)

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
      do IK=1,4
        COB(1,IK) = CK(IK,4)*CA3
        COB(2,IK) = CK(IK,4)*CC3
        COB(3,IK) = CK(IK,4)*CB3
        COB(4,IK) = CK(IK,2)*CA1+CK(IK,3)*CA2
        COB(5,IK) = CK(IK,2)*CB1+CK(IK,3)*CB2
        COB(6,IK) = CK(IK,1)
      end do

      ! Function values
      do i=1,npointsfunc   ! either 0 or npoints

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
      end do

      ! x-derivatives

      do i=1,npointsderx   ! either 0 or npoints
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
      end do

      ! y-derivatives

      do i=1,npointsdery   ! either 0 or npoints
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

      end do

    end do ! ielement
    !$omp end parallel do

  end if

  contains

    ! Auxiliary functions

    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function

    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function

    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function

    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(in) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine

!**************************************************************************
! Element subroutines for parametric Q1~ element, midpoint value based.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>

  pure subroutine elem_E031 (celement, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D) :: Dhelp
  real(DP) :: dx,dy

  real(DP) :: dxj !auxiliary variable

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
  ! That is even faster than when using three IF commands for preventing
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

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_E031_mult (celement, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC,i) = 0.25E0_DP*&
          (-(Dpoints(1,i)**2-Dpoints(2,i)**2)-2.0_DP*Dpoints(2,i)+1.0_DP)
      Dbas(2,DER_FUNC,i) = 0.25E0_DP*&
          ( (Dpoints(1,i)**2-Dpoints(2,i)**2)+2.0_DP*Dpoints(1,i)+1.0_DP)
      Dbas(3,DER_FUNC,i) = 0.25E0_DP*&
          (-(Dpoints(1,i)**2-Dpoints(2,i)**2)+2.0_DP*Dpoints(2,i)+1.0_DP)
      Dbas(4,DER_FUNC,i) = 0.25E0_DP*&
          ( (Dpoints(1,i)**2-Dpoints(2,i)**2)-2.0_DP*Dpoints(1,i)+1.0_DP)
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj(:) = 0.5E0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
      Dhelp(1,1,i) = -Dpoints(1,i)
      Dhelp(2,1,i) =  Dpoints(1,i)+1.0_DP
      Dhelp(3,1,i) = -Dpoints(1,i)
      Dhelp(4,1,i) =  Dpoints(1,i)-1.0_DP
      Dhelp(1,2,i) =  Dpoints(2,i)-1.0_DP
      Dhelp(2,2,i) = -Dpoints(2,i)
      Dhelp(3,2,i) =  Dpoints(2,i)+1.0_DP
      Dhelp(4,2,i) = -Dpoints(2,i)
    end do

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_E031_sim (celement, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, &
                            Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 0.25E0_DP*&
            (-(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-2.0_DP*Dpoints(2,i,j)+1.0_DP)
        Dbas(2,DER_FUNC,i,j) = 0.25E0_DP*&
            ( (Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+2.0_DP*Dpoints(1,i,j)+1.0_DP)
        Dbas(3,DER_FUNC,i,j) = 0.25E0_DP*&
            (-(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+2.0_DP*Dpoints(2,i,j)+1.0_DP)
        Dbas(4,DER_FUNC,i,j) = 0.25E0_DP*&
            ( (Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-2.0_DP*Dpoints(1,i,j)+1.0_DP)
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj,Dhelp) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 0.5E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
        Dhelp(1,1,i) = -Dpoints(1,i,j)
        Dhelp(2,1,i) =  Dpoints(1,i,j)+1.0_DP
        Dhelp(3,1,i) = -Dpoints(1,i,j)
        Dhelp(4,1,i) =  Dpoints(1,i,j)-1.0_DP
        Dhelp(1,2,i) =  Dpoints(2,i,j)-1.0_DP
        Dhelp(2,2,i) = -Dpoints(2,i,j)
        Dhelp(3,2,i) =  Dpoints(2,i,j)+1.0_DP
        Dhelp(4,2,i) = -Dpoints(2,i,j)
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(1,1,i) &
                                    - Djac(2,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(2,1,i) &
                                    - Djac(2,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(3,1,i) &
                                    - Djac(2,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(4,1,i) &
                                    - Djac(2,i,j) * Dhelp(4,2,i))
        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(1,1,i) &
                                    - Djac(1,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(2,1,i) &
                                    - Djac(1,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(3,1,i) &
                                    - Djac(1,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(4,1,i) &
                                    - Djac(1,i,j) * Dhelp(4,2,i))
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric Q2~ element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>

  pure subroutine elem_E050 (celement, Dcoords, itwistIndex, Djac, ddetj, &
                             Bder, Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2T.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges.
  integer(I32), intent(in) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(9,2) :: Dhelp

  ! Scaling factors for basis functions 5..8
  real(DP) :: d5,d6,d7,d8

  ! auxiliary variables
  real(DP) :: dx,dy,dxj

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

  ! The Q2~ element is specified by nine polynomials on the reference element.
  ! These nine polynomials are:
  !
  !  p_1(x,y) =  3/4*y*( x^2 + y - 1) - 1/4
  !  p_2(x,y) =  3/4*x*(-y^2 + x + 1) - 1/4
  !  p_3(x,y) =  3/4*y*(-x^2 + y + 1) - 1/4
  !  p_4(x,y) =  3/4*x*( y^2 + x - 1) - 1/4
  !  p_5(x,y) = -3/8*x*(y*(y*( 5*y - 6) - 5*x^2 + 2) + 2)
  !  p_6(x,y) = -3/8*y*(x*(x*(-5*x - 6) + 5*y^2 - 2) + 2)
  !  p_7(x,y) = -3/8*x*(y*(y*( 5*y + 6) - 5*x^2 + 2) - 2)
  !  p_8(x,y) = -3/8*y*(x*(x*(-5*x + 6) + 5*y^2 - 2) - 2)
  !  p_9(x,y) = -3/2*(x^2 + y^2) + 2
  !
  ! Remark: Since the degree of the monoms goes up to 3 for these basis
  ! polynomials, the polynomials are evaluated using the Horner scheme.

  ! Clear the output array
  !Dbas = 0.0_DP
  dx = Dpoint(1)
  dy = Dpoint(2)

  ! Get the twist indices that define the orientation of our edges.
  ! A value of 1 is standard, a value of -1 results in a change of the
  ! sign in the basis functions p_5, ..., p_8.
  ! itwistIndex is a bitfield. Each bit specifies the orientation of an edge.
  ! We use the bit to calculate a "1.0" if the edge has positive orientation
  ! and "-1.0" if it has negative orientation.
  d5 = real(1-iand(int(ishft(itwistIndex, 1)),2),DP)
  d6 = real(1-iand(int(ishft(itwistIndex, 0)),2),DP)
  d7 = real(1-iand(int(ishft(itwistIndex,-1)),2),DP)
  d8 = real(1-iand(int(ishft(itwistIndex,-2)),2),DP)

  ! Remark: The Q2~-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) =  Q3*dy*( dx**2 + dy - P1) - Q1
    Dbas(2,DER_FUNC) =  Q3*dx*(-dy**2 + dx + P1) - Q1
    Dbas(3,DER_FUNC) =  Q3*dy*(-dx**2 + dy + P1) - Q1
    Dbas(4,DER_FUNC) =  Q3*dx*( dy**2 + dx - P1) - Q1
    Dbas(5,DER_FUNC) = -Q2*dx*(dy*(dy*( P4*dy - P5) - P4*dx**2 + P2) + P2)*d5
    Dbas(6,DER_FUNC) = -Q2*dy*(dx*(dx*(-P4*dx - P5) + P4*dy**2 - P2) + P2)*d6
    Dbas(7,DER_FUNC) = -Q2*dx*(dy*(dy*( P4*dy + P5) - P4*dx**2 + P2) - P2)*d7
    Dbas(8,DER_FUNC) = -Q2*dy*(dx*(dx*(-P4*dx + P5) + P4*dy**2 - P2) - P2)*d8
    Dbas(9,DER_FUNC) = -Q4*(dx**2 + dy**2) + P2

    ! old FEAT 1.3 code
!    Dbas(1,DER_FUNC) = -1.D0/4.D0-y/2+3.D0/4.D0*y**2+(-1.D0/2.D0+3.D0/2.D0*x**2)*y/2
!    Dbas(2,DER_FUNC) = -1.D0/4.D0+x/2+3.D0/4.D0*x**2-(-1.D0/2.D0+3.D0/2.D0*y**2)*x/2
!    Dbas(3,DER_FUNC) = -1.D0/4.D0+y/2+3.D0/4.D0*y**2-(-1.D0/2.D0+3.D0/2.D0*x**2)*y/2
!    Dbas(4,DER_FUNC) = -1.D0/4.D0-x/2+3.D0/4.D0*x**2+(-1.D0/2.D0+3.D0/2.D0*y**2)*x/2
!    Dbas(5,DER_FUNC) = (-3.D0/4.D0*x*y+3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*y**2)*x+3.D0/&
!      4.D0*(5.D0/2.D0*x**3-3.D0/2.D0*x)*y-3.D0/4.D0*x*(5.D0/2.D0*y**3-3.D0/2.D0*y))*d5
!    Dbas(6,DER_FUNC) = (3.D0/4.D0*x*y+3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*x**2)*y+3.D0/&
!      4.D0*(5.D0/2.D0*x**3-3.D0/2.D0*x)*y-3.D0/4.D0*x*(5.D0/2.D0*y**3-3.D0/2.D0*y))*d6
!    Dbas(7,DER_FUNC) = (-3.D0/4.D0*x*y-3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*y**2)*x+3.D0/&
!      4.D0*(5.D0/2.D0*x**3-3.D0/2.D0*x)*y-3.D0/4.D0*x*(5.D0/2.D0*y**3-3.D0/2.D0*y))*d7
!    Dbas(8,DER_FUNC) = (3.D0/4.D0*x*y-3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*x**2)*y+3.D0/&
!      4.D0*(5.D0/2.D0*x**3-3.D0/2.D0*x)*y-3.D0/4.D0*x*(5.D0/2.D0*y**3-3.D0/2.D0*y))*d8
!    Dbas(9,DER_FUNC) = 2.D0-3.D0/2.D0*x**2-3.D0/2.D0*y**2


!  endif

  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 1.0_DP / ddetj

    ! x- and y-derivatives on reference element
    Dhelp(1,1) =  Q4*dx*dy
    Dhelp(2,1) = -Q3*dy**2 + Q4*dx + Q3
    Dhelp(3,1) = -Q4*dx*dy
    Dhelp(4,1) =  Q3*dy**2 + Q4*dx - Q3
    Dhelp(5,1) =  dy*(dy*(-Q6*dy + Q5) + Q7*dx**2 - Q3) - Q3
    Dhelp(6,1) = -Q2*dy*(dx*(-P7*dx - P6) + P4*dy**2 - P2)
    Dhelp(7,1) =  dy*(dy*(-Q6*dy - Q5) + Q7*dx**2 - Q3) + Q3
    Dhelp(8,1) = -Q2*dy*(dx*(-P7*dx + P6) + P4*dy**2 - P2)
    Dhelp(9,1) = -P3*dx
    Dhelp(1,2) =  Q3*dx**2 + Q4*dy - Q3
    Dhelp(2,2) = -Q4*dx*dy
    Dhelp(3,2) = -Q3*dx**2 + Q4*dy + Q3
    Dhelp(4,2) =  Q4*dx*dy
    Dhelp(5,2) = -Q2*dx*(dy*( P7*dy - P6) - P4*dx**2 + P2)
    Dhelp(6,2) =  dx*(dx*( Q6*dx + Q5) - Q7*dy**2 + Q3) - Q3
    Dhelp(7,2) =  Q2*dx*(dy*(-P7*dy - P6) + P4*dx**2 - P2)
    Dhelp(8,2) =  dx*(dx*( Q6*dx - Q5) - Q7*dy**2 + Q3) + Q3
    Dhelp(9,2) = -P3*dy
    ! 'old' FEAT 1.3 code
!    Dhelp(1,1) = 3.D0/2.D0*x*y
!    Dhelp(2,1) = 3.D0/4.D0+3.D0/2.D0*x-3.D0/4.D0*y**2
!    Dhelp(3,1) = -3.D0/2.D0*x*y
!    Dhelp(4,1) = -3.D0/4.D0+3.D0/2.D0*x+3.D0/4.D0*y**2
!    Dhelp(5,1) = 3.D0/8.D0*y-3.D0/4.D0+9.D0/4.D0*y**2+3.D0/4.D0*(15.D0/2.D0*&
!                 x**2-3.D0/2.D0)*y-15.D0/8.D0*y**3
!    Dhelp(6,1) = 15.D0/8.D0*y+9.D0/2.D0*x*y+3.D0/4.D0*(15.D0/2.D0*x**2-3.D0/&
!                 2.D0)*y-15.D0/8.D0*y**3
!    Dhelp(7,1) = 3.D0/8.D0*y+3.D0/4.D0-9.D0/4.D0*y**2+3.D0/4.D0*(15.D0/2.D0*&
!                 x**2-3.D0/2.D0)*y-15.D0/8.D0*y**3
!    Dhelp(8,1) = 15.D0/8.D0*y-9.D0/2.D0*x*y+3.D0/4.D0*(15.D0/2.D0*x**2-3.D0/&
!                 2.D0)*y-15.D0/8.D0*y**3
!    Dhelp(9,1) = -3.D0*x
!    Dhelp(1,2) = -3.D0/4.D0+3.D0/2.D0*y+3.D0/4.D0*x**2
!    Dhelp(2,2) = -3.D0/2.D0*x*y
!    Dhelp(3,2) = 3.D0/4.D0+3.D0/2.D0*y-3.D0/4.D0*x**2
!    Dhelp(4,2) = 3.D0/2.D0*x*y
!    Dhelp(5,2) = -15.D0/8.D0*x+9.D0/2.D0*x*y+15.D0/8.D0*x**3-3.D0/4.D0*x*&
!                 (15.D0/2.D0*y**2-3.D0/2.D0)
!    Dhelp(6,2) = -3.D0/8.D0*x-3.D0/4.D0+9.D0/4.D0*x**2+15.D0/8.D0*x**3-3.D0/&
!                 4.D0*x*(15.D0/2.D0*y**2-3.D0/2.D0)
!    Dhelp(7,2) = -15.D0/8.D0*x-9.D0/2.D0*x*y+15.D0/8.D0*x**3-3.D0/4.D0*x*&
!                 (15.D0/2.D0*y**2-3.D0/2.D0)
!    Dhelp(8,2) = -3.D0/8.D0*x+3.D0/4.D0-9.D0/4.D0*x**2+15.D0/8.D0*x**3-3.D0/&
!                 4.D0*x*(15.D0/2.D0*y**2-3.D0/2.D0)
!    Dhelp(9,2) = -3*y

    ! x-derivatives on current element
!    if (Bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(1,1) - Djac(2)*Dhelp(1,2))
      Dbas(2,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(2,1) - Djac(2)*Dhelp(2,2))
      Dbas(3,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(3,1) - Djac(2)*Dhelp(3,2))
      Dbas(4,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(4,1) - Djac(2)*Dhelp(4,2))
      Dbas(5,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(5,1) - Djac(2)*Dhelp(5,2))*d5
      Dbas(6,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(6,1) - Djac(2)*Dhelp(6,2))*d6
      Dbas(7,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(7,1) - Djac(2)*Dhelp(7,2))*d7
      Dbas(8,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(8,1) - Djac(2)*Dhelp(8,2))*d8
      Dbas(9,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(9,1) - Djac(2)*Dhelp(9,2))
!    endif

    ! y-derivatives on current element
!    if (Bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(1,1) - Djac(1)*Dhelp(1,2))
      Dbas(2,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(2,1) - Djac(1)*Dhelp(2,2))
      Dbas(3,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(3,1) - Djac(1)*Dhelp(3,2))
      Dbas(4,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(4,1) - Djac(1)*Dhelp(4,2))
      Dbas(5,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(5,1) - Djac(1)*Dhelp(5,2))*d5
      Dbas(6,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(6,1) - Djac(1)*Dhelp(6,2))*d6
      Dbas(7,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(7,1) - Djac(1)*Dhelp(7,2))*d7
      Dbas(8,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(8,1) - Djac(1)*Dhelp(8,2))*d8
      Dbas(9,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(9,1) - Djac(1)*Dhelp(9,2))
!    endif
!  endif

  end subroutine


  !************************************************************************

!<subroutine>

  pure subroutine elem_E050_mult (celement, Dcoords, itwistIndex, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2T.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges.
  integer(I32), intent(in) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(9,2,npoints) :: Dhelp

  ! Scaling factors for basis functions 5..8
  real(DP) :: d5,d6,d7,d8

  ! Auxiliary variables
  real(DP), dimension(npoints) :: dxj
  real(DP) :: dx,dy

  integer :: i   ! point counter

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Get the twist indices for this element
  ! itwistIndex is a bitfield. Each bit specifies the orientation of an edge.
  ! We use the bit to calculate a "1.0" if the edge has positive orientation
  ! and "-1.0" if it has negative orientation.
  d5 = real(1-iand(int(ishft(itwistIndex, 1)),2),DP)
  d6 = real(1-iand(int(ishft(itwistIndex, 0)),2),DP)
  d7 = real(1-iand(int(ishft(itwistIndex,-1)),2),DP)
  d8 = real(1-iand(int(ishft(itwistIndex,-2)),2),DP)

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dbas(1,DER_FUNC,i) =  Q3*dy*( dx**2 + dy - P1) - Q1
      Dbas(2,DER_FUNC,i) =  Q3*dx*(-dy**2 + dx + P1) - Q1
      Dbas(3,DER_FUNC,i) =  Q3*dy*(-dx**2 + dy + P1) - Q1
      Dbas(4,DER_FUNC,i) =  Q3*dx*( dy**2 + dx - P1) - Q1
      Dbas(5,DER_FUNC,i) = -Q2*dx*(dy*(dy*( P4*dy - P5) - P4*dx**2 + P2) + P2)*d5
      Dbas(6,DER_FUNC,i) = -Q2*dy*(dx*(dx*(-P4*dx - P5) + P4*dy**2 - P2) + P2)*d6
      Dbas(7,DER_FUNC,i) = -Q2*dx*(dy*(dy*( P4*dy + P5) - P4*dx**2 + P2) - P2)*d7
      Dbas(8,DER_FUNC,i) = -Q2*dy*(dx*(dx*(-P4*dx + P5) + P4*dy**2 - P2) - P2)*d8
      Dbas(9,DER_FUNC,i) = -Q4*(dx**2 + dy**2) + P2
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN

    !x- and y-derivatives on reference element
    do i=1,npoints
      dxj(i) = 1.0_DP / Ddetj(i)
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dhelp(1,1,i) =  Q4*dx*dy
      Dhelp(2,1,i) = -Q3*dy**2 + Q4*dx + Q3
      Dhelp(3,1,i) = -Q4*dx*dy
      Dhelp(4,1,i) =  Q3*dy**2 + Q4*dx - Q3
      Dhelp(5,1,i) =  dy*(dy*(-Q6*dy + Q5) + Q7*dx**2 - Q3) - Q3
      Dhelp(6,1,i) = -Q2*dy*(dx*(-P7*dx - P6) + P4*dy**2 - P2)
      Dhelp(7,1,i) =  dy*(dy*(-Q6*dy - Q5) + Q7*dx**2 - Q3) + Q3
      Dhelp(8,1,i) = -Q2*dy*(dx*(-P7*dx + P6) + P4*dy**2 - P2)
      Dhelp(9,1,i) = -P3*dx
      Dhelp(1,2,i) =  Q3*dx**2 + Q4*dy - Q3
      Dhelp(2,2,i) = -Q4*dx*dy
      Dhelp(3,2,i) = -Q3*dx**2 + Q4*dy + Q3
      Dhelp(4,2,i) =  Q4*dx*dy
      Dhelp(5,2,i) = -Q2*dx*(dy*( P7*dy - P6) - P4*dx**2 + P2)
      Dhelp(6,2,i) =  dx*(dx*( Q6*dx + Q5) - Q7*dy**2 + Q3) - Q3
      Dhelp(7,2,i) =  Q2*dx*(dy*(-P7*dy - P6) + P4*dx**2 - P2)
      Dhelp(8,2,i) =  dx*(dx*( Q6*dx - Q5) - Q7*dy**2 + Q3) + Q3
      Dhelp(9,2,i) = -P3*dy
    end do

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        Dbas(1,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(1,1,i) &
                                         - Djac(2,i)*Dhelp(1,2,i))
        Dbas(2,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(2,1,i) &
                                         - Djac(2,i)*Dhelp(2,2,i))
        Dbas(3,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(3,1,i) &
                                         - Djac(2,i)*Dhelp(3,2,i))
        Dbas(4,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(4,1,i) &
                                         - Djac(2,i)*Dhelp(4,2,i))
        Dbas(5,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(5,1,i) &
                                         - Djac(2,i)*Dhelp(5,2,i))*d5
        Dbas(6,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(6,1,i) &
                                         - Djac(2,i)*Dhelp(6,2,i))*d6
        Dbas(7,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(7,1,i) &
                                         - Djac(2,i)*Dhelp(7,2,i))*d7
        Dbas(8,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(8,1,i) &
                                         - Djac(2,i)*Dhelp(8,2,i))*d8
        Dbas(9,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(9,1,i) &
                                         - Djac(2,i)*Dhelp(9,2,i))
!      END DO
!    ENDIF

    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(1,1,i) &
                                         - Djac(1,i)*Dhelp(1,2,i))
        Dbas(2,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(2,1,i) &
                                         - Djac(1,i)*Dhelp(2,2,i))
        Dbas(3,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(3,1,i) &
                                         - Djac(1,i)*Dhelp(3,2,i))
        Dbas(4,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(4,1,i) &
                                         - Djac(1,i)*Dhelp(4,2,i))
        Dbas(5,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(5,1,i) &
                                         - Djac(1,i)*Dhelp(5,2,i))*d5
        Dbas(6,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(6,1,i) &
                                         - Djac(1,i)*Dhelp(6,2,i))*d6
        Dbas(7,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(7,1,i) &
                                         - Djac(1,i)*Dhelp(7,2,i))*d7
        Dbas(8,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(8,1,i) &
                                         - Djac(1,i)*Dhelp(8,2,i))*d8
        Dbas(9,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(9,1,i) &
                                         - Djac(1,i)*Dhelp(9,2,i))
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_E050_sim (celement, Dcoords, ItwistIndex, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2T.
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

  ! List of twist indices. For every edge/face on every cell, the twist
  ! index defines the orientation of the edge/face.
  ! Array with DIMENSION(1:NVE/NVA,nelements)
  integer(I32), dimension(:), intent(in) :: ItwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(9,2,npoints) :: Dhelp

  ! Scaling factors for basis functions 5..8
  real(DP) :: d5,d6,d7,d8

  ! Auxiliary variables
  real(DP), dimension(npoints) :: dxj
  real(DP) :: dx,dy

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i,dx,dy,d5,d6,d7,d8) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      ! Get the twist indices for this element.
      ! ItwistIndex(.) is a bitfield. Each bit specifies the orientation of an edge.
      ! We use the bit to calculate a "1.0" if the edge has positive orientation
      ! and "-1.0" if it has negative orientation.
      d5 = real(1-iand(int(ishft(ItwistIndex(j), 1)),2),DP)
      d6 = real(1-iand(int(ishft(ItwistIndex(j), 0)),2),DP)
      d7 = real(1-iand(int(ishft(ItwistIndex(j),-1)),2),DP)
      d8 = real(1-iand(int(ishft(ItwistIndex(j),-2)),2),DP)

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dbas(1,DER_FUNC,i,j) =  Q3*dy*( dx**2 + dy - P1) - Q1
        Dbas(2,DER_FUNC,i,j) =  Q3*dx*(-dy**2 + dx + P1) - Q1
        Dbas(3,DER_FUNC,i,j) =  Q3*dy*(-dx**2 + dy + P1) - Q1
        Dbas(4,DER_FUNC,i,j) =  Q3*dx*( dy**2 + dx - P1) - Q1
        Dbas(5,DER_FUNC,i,j) = -Q2*dx*(dy*(dy*( P4*dy - P5) - P4*dx**2 + P2) + P2)*d5
        Dbas(6,DER_FUNC,i,j) = -Q2*dy*(dx*(dx*(-P4*dx - P5) + P4*dy**2 - P2) + P2)*d6
        Dbas(7,DER_FUNC,i,j) = -Q2*dx*(dy*(dy*( P4*dy + P5) - P4*dx**2 + P2) - P2)*d7
        Dbas(8,DER_FUNC,i,j) = -Q2*dy*(dx*(dx*(-P4*dx + P5) + P4*dy**2 - P2) - P2)*d8
        Dbas(9,DER_FUNC,i,j) = -Q4*(dx**2 + dy**2) + P2
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj,Dhelp,dx,dy,d5,d6,d7,d8) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      ! Get the twist indices for this element
      ! ItwistIndex(.) is a bitfield. Each bit specifies the orientation of an edge.
      ! We use the bit to calculate a "1.0" if the edge has positive orientation
      ! and "-1.0" if it has negative orientation.
      d5 = real(1-iand(int(ishft(ItwistIndex(j), 1)),2),DP)
      d6 = real(1-iand(int(ishft(ItwistIndex(j), 0)),2),DP)
      d7 = real(1-iand(int(ishft(ItwistIndex(j),-1)),2),DP)
      d8 = real(1-iand(int(ishft(ItwistIndex(j),-2)),2),DP)

      !x- and y-derivatives on reference element
      do i=1,npoints
        dxj(i) = 1.0_DP / Ddetj(i,j)
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dhelp(1,1,i) =  Q4*dx*dy
        Dhelp(2,1,i) = -Q3*dy**2 + Q4*dx + Q3
        Dhelp(3,1,i) = -Q4*dx*dy
        Dhelp(4,1,i) =  Q3*dy**2 + Q4*dx - Q3
        Dhelp(5,1,i) =  dy*(dy*(-Q6*dy + Q5) + Q7*dx**2 - Q3) - Q3
        Dhelp(6,1,i) = -Q2*dy*(dx*(-P7*dx - P6) + P4*dy**2 - P2)
        Dhelp(7,1,i) =  dy*(dy*(-Q6*dy - Q5) + Q7*dx**2 - Q3) + Q3
        Dhelp(8,1,i) = -Q2*dy*(dx*(-P7*dx + P6) + P4*dy**2 - P2)
        Dhelp(9,1,i) = -P3*dx
        Dhelp(1,2,i) =  Q3*dx**2 + Q4*dy - Q3
        Dhelp(2,2,i) = -Q4*dx*dy
        Dhelp(3,2,i) = -Q3*dx**2 + Q4*dy + Q3
        Dhelp(4,2,i) =  Q4*dx*dy
        Dhelp(5,2,i) = -Q2*dx*(dy*( P7*dy - P6) - P4*dx**2 + P2)
        Dhelp(6,2,i) =  dx*(dx*( Q6*dx + Q5) - Q7*dy**2 + Q3) - Q3
        Dhelp(7,2,i) =  Q2*dx*(dy*(-P7*dy - P6) + P4*dx**2 - P2)
        Dhelp(8,2,i) =  dx*(dx*( Q6*dx - Q5) - Q7*dy**2 + Q3) + Q3
        Dhelp(9,2,i) = -P3*dy
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(1,1,i) &
                                           - Djac(2,i,j)*Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(2,1,i) &
                                           - Djac(2,i,j)*Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(3,1,i) &
                                           - Djac(2,i,j)*Dhelp(3,2,i))
          Dbas(4,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(4,1,i) &
                                           - Djac(2,i,j)*Dhelp(4,2,i))
          Dbas(5,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(5,1,i) &
                                           - Djac(2,i,j)*Dhelp(5,2,i))*d5
          Dbas(6,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(6,1,i) &
                                           - Djac(2,i,j)*Dhelp(6,2,i))*d6
          Dbas(7,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(7,1,i) &
                                           - Djac(2,i,j)*Dhelp(7,2,i))*d7
          Dbas(8,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(8,1,i) &
                                           - Djac(2,i,j)*Dhelp(8,2,i))*d8
          Dbas(9,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(9,1,i) &
                                           - Djac(2,i,j)*Dhelp(9,2,i))
        !end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        !do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(1,1,i) &
                                           - Djac(1,i,j)*Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(2,1,i) &
                                           - Djac(1,i,j)*Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(3,1,i) &
                                           - Djac(1,i,j)*Dhelp(3,2,i))
          Dbas(4,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(4,1,i) &
                                           - Djac(1,i,j)*Dhelp(4,2,i))
          Dbas(5,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(5,1,i) &
                                           - Djac(1,i,j)*Dhelp(5,2,i))*d5
          Dbas(6,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(6,1,i) &
                                           - Djac(1,i,j)*Dhelp(6,2,i))*d6
          Dbas(7,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(7,1,i) &
                                           - Djac(1,i,j)*Dhelp(7,2,i))*d7
          Dbas(8,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(8,1,i) &
                                           - Djac(1,i,j)*Dhelp(8,2,i))*d8
          Dbas(9,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(9,1,i) &
                                           - Djac(1,i,j)*Dhelp(9,2,i))
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric Q2~ element with bubble.
!
! The element is nearly the same as Q2~, but adds an additional bubble
! to the basis functions. This should archive the element to be stable
! against mesh distortion.
!
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>

  pure subroutine elem_EB50 (celement, Dcoords, itwistIndex, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2TB.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges.
  integer(I32), intent(in) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(10,2) :: Dhelp

  ! Scaling factors for basis functions 5..8
  real(DP) :: d5,d6,d7,d8

  ! auxiliary variables
  real(DP) :: dx,dy,dx2,dy2,dxj

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: Q8 = 6.25_DP  ! = 25/4
  real(DP), parameter :: Q9 = 37.5_DP  ! = 75/2
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

  ! The Q2~ element with bubble is specified by ten polynomials on the
  ! reference element. These ten polynomials are:
  !
  !  p_1 (x,y) =  3/4*y*( x^2 + y - 1) - 1/4
  !  p_2 (x,y) =  3/4*x*(-y^2 + x + 1) - 1/4
  !  p_3 (x,y) =  3/4*y*(-x^2 + y + 1) - 1/4
  !  p_4 (x,y) =  3/4*x*( y^2 + x - 1) - 1/4
  !  p_5 (x,y) = -3/8*x*(y*(y*( 5*y - 6) - 5*x^2 + 2) + 2)
  !  p_6 (x,y) = -3/8*y*(x*(x*(-5*x - 6) + 5*y^2 - 2) + 2)
  !  p_7 (x,y) = -3/8*x*(y*(y*( 5*y + 6) - 5*x^2 + 2) - 2)
  !  p_8 (x,y) = -3/8*y*(x*(x*(-5*x + 6) + 5*y^2 - 2) - 2)
  !  p_9 (x,y) = -3/2*(x^2 + y^2) + 2
  !  p_10(x,y) = 25/4*(3*x^2 - 1)*(3*y^2 - 1)
  !
  ! These ten polynomials are constructed in a way such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,10
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (p_i(Ej(t))      ) d(t) = kronecker(i,j  ) * |ej|
  !     Int_[-1,1] (p_i(Ej(t))*L1(t)) d(t) = kronecker(i,j+4) * |ej|
  !   }
  !   Int_T (p_i(x,y)            ) d(x,y) = kronecker(i, 9) * |T|
  !   Int_T (p_i(x,y)*L2(x)*L2(y)) d(x,y) = kronecker(i,10) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th edge of the reference quadrilateral
  ! Ej: [-1,1] -> ej being the parametrisation of an edge ej
  ! |ej| = 2 being the length of the edge ej
  ! T being the reference quadrilateral
  ! |T| = 4 being the area of the reference quadrilateral
  ! L1 and L2 being the first two Legendre-Polynomials:
  ! L1(x) := x
  ! L2(x) := 1/2*(3*x^2 - 1)
  !
  ! Remark: Since the degree of the monoms goes up to 3 for these basis
  ! polynomials, the polynomials are evaluated using the Horner scheme.

  ! Clear the output array
  !Dbas = 0.0_DP
  dx = Dpoint(1)
  dy = Dpoint(2)
  dx2 = dx**2
  dy2 = dy**2

  ! Get the twist indices that define the orientation of our edges.
  ! A value of 1 is standard, a value of -1 results in a change of the
  ! sign in the basis functions p_5, ..., p_8.
  ! itwistIndex(.) is a bitfield. Each bit specifies the orientation of an edge.
  ! We use the bit to calculate a "1.0" if the edge has positive orientation
  ! and "-1.0" if it has negative orientation.
  d5 = real(1-iand(int(ishft(itwistIndex, 1)),2),DP)
  d6 = real(1-iand(int(ishft(itwistIndex, 0)),2),DP)
  d7 = real(1-iand(int(ishft(itwistIndex,-1)),2),DP)
  d8 = real(1-iand(int(ishft(itwistIndex,-2)),2),DP)

  ! Remark: The Q2~-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) =  Q3*dy*( dx2 + dy - P1) - Q1
    Dbas(2,DER_FUNC) =  Q3*dx*(-dy2 + dx + P1) - Q1
    Dbas(3,DER_FUNC) =  Q3*dy*(-dx2 + dy + P1) - Q1
    Dbas(4,DER_FUNC) =  Q3*dx*( dy2 + dx - P1) - Q1
    Dbas(5,DER_FUNC) = -Q2*dx*(dy*(dy*( P4*dy - P5) - P4*dx2 + P2) + P2)*d5
    Dbas(6,DER_FUNC) = -Q2*dy*(dx*(dx*(-P4*dx - P5) + P4*dy2 - P2) + P2)*d6
    Dbas(7,DER_FUNC) = -Q2*dx*(dy*(dy*( P4*dy + P5) - P4*dx2 + P2) - P2)*d7
    Dbas(8,DER_FUNC) = -Q2*dy*(dx*(dx*(-P4*dx + P5) + P4*dy2 - P2) - P2)*d8
    Dbas(9,DER_FUNC) = -Q4*(dx2 + dy2) + P2
    Dbas(10,DER_FUNC)=  Q8*(P3*dx2 - P1)*(P3*dy2 - P1)

    ! old FEAT 1.3 code
!    Dbas(1,DER_FUNC) = -1.D0/4.D0-dy/2+3.D0/4.D0*dy**2+(-1.D0/2.D0+3.D0/2.D0*dx**2)*dy/2
!    Dbas(2,DER_FUNC) = -1.D0/4.D0+dx/2+3.D0/4.D0*dx**2-(-1.D0/2.D0+3.D0/2.D0*dy**2)*dx/2
!    Dbas(3,DER_FUNC) = -1.D0/4.D0+dy/2+3.D0/4.D0*dy**2-(-1.D0/2.D0+3.D0/2.D0*dx**2)*dy/2
!    Dbas(4,DER_FUNC) = -1.D0/4.D0-dx/2+3.D0/4.D0*dx**2+(-1.D0/2.D0+3.D0/2.D0*dy**2)*dx/2
!    Dbas(5,DER_FUNC) = (-3.D0/4.D0*dx*dy+3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*dy**2)*dx+3.D0/&
!      4.D0*(5.D0/2.D0*dx**3-3.D0/2.D0*dx)*dy-3.D0/4.D0*dx*(5.D0/2.D0*dy**3-3.D0/2.D0*dy))*d5
!    Dbas(6,DER_FUNC) = (3.D0/4.D0*dx*dy+3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*dx**2)*dy+3.D0/&
!      4.D0*(5.D0/2.D0*dx**3-3.D0/2.D0*dx)*dy-3.D0/4.D0*dx*(5.D0/2.D0*dy**3-3.D0/2.D0*dy))*d6
!    Dbas(7,DER_FUNC) = (-3.D0/4.D0*dx*dy-3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*dy**2)*dx+3.D0/&
!      4.D0*(5.D0/2.D0*dx**3-3.D0/2.D0*dx)*dy-3.D0/4.D0*dx*(5.D0/2.D0*dy**3-3.D0/2.D0*dy))*d7
!    Dbas(8,DER_FUNC) = (3.D0/4.D0*dx*dy-3.D0/2.D0*(-1.D0/2.D0+3.D0/2.D0*dx**2)*dy+3.D0/&
!      4.D0*(5.D0/2.D0*dx**3-3.D0/2.D0*dx)*dy-3.D0/4.D0*dx*(5.D0/2.D0*dy**3-3.D0/2.D0*dy))*d8
!    Dbas(9,DER_FUNC) = 2.D0-3.D0/2.D0*dx**2-3.D0/2.D0*dy**2
!    Dbas(10,DER_FUNC) = (-1.D0/2.D0+3.D0/2.D0*dx**2)*(-1.D0/2.D0+3.D0/2.D0*dy**2)


!  endif

  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 1.0_DP / ddetj

    ! x- and y-derivatives on reference element
    Dhelp(1,1) =  Q4*dx*dy
    Dhelp(2,1) = -Q3*dy2 + Q4*dx + Q3
    Dhelp(3,1) = -Q4*dx*dy
    Dhelp(4,1) =  Q3*dy2 + Q4*dx - Q3
    Dhelp(5,1) =  dy*(dy*(-Q6*dy + Q5) + Q7*dx2 - Q3) - Q3
    Dhelp(6,1) = -Q2*dy*(dx*(-P7*dx - P6) + P4*dy2 - P2)
    Dhelp(7,1) =  dy*(dy*(-Q6*dy - Q5) + Q7*dx2 - Q3) + Q3
    Dhelp(8,1) = -Q2*dy*(dx*(-P7*dx + P6) + P4*dy2 - P2)
    Dhelp(9,1) = -P3*dx
    Dhelp(10,1)=  Q9*dx*(P3*dy2 - P1)
    Dhelp(1,2) =  Q3*dx2 + Q4*dy - Q3
    Dhelp(2,2) = -Q4*dx*dy
    Dhelp(3,2) = -Q3*dx2 + Q4*dy + Q3
    Dhelp(4,2) =  Q4*dx*dy
    Dhelp(5,2) = -Q2*dx*(dy*( P7*dy - P6) - P4*dx2 + P2)
    Dhelp(6,2) =  dx*(dx*( Q6*dx + Q5) - Q7*dy2 + Q3) - Q3
    Dhelp(7,2) =  Q2*dx*(dy*(-P7*dy - P6) + P4*dx2 - P2)
    Dhelp(8,2) =  dx*(dx*( Q6*dx - Q5) - Q7*dy2 + Q3) + Q3
    Dhelp(9,2) = -P3*dy
    Dhelp(10,2)=  Q9*dy*(P3*dx2 - P1)

    ! 'old' FEAT 1.3 code
!    Dhelp(1,1) = 3.D0/2.D0*dx*dy
!    Dhelp(2,1) = 3.D0/4.D0+3.D0/2.D0*dx-3.D0/4.D0*dy**2
!    Dhelp(3,1) = -3.D0/2.D0*dx*dy
!    Dhelp(4,1) = -3.D0/4.D0+3.D0/2.D0*dx+3.D0/4.D0*dy**2
!    Dhelp(5,1) = 3.D0/8.D0*dy-3.D0/4.D0+9.D0/4.D0*dy**2+3.D0/4.D0*(15.D0/2.D0*&
!                 dx**2-3.D0/2.D0)*dy-15.D0/8.D0*dy**3
!    Dhelp(6,1) = 15.D0/8.D0*dy+9.D0/2.D0*dx*dy+3.D0/4.D0*(15.D0/2.D0*dx**2-3.D0/&
!                 2.D0)*dy-15.D0/8.D0*dy**3
!    Dhelp(7,1) = 3.D0/8.D0*dy+3.D0/4.D0-9.D0/4.D0*dy**2+3.D0/4.D0*(15.D0/2.D0*&
!                 dx**2-3.D0/2.D0)*dy-15.D0/8.D0*dy**3
!    Dhelp(8,1) = 15.D0/8.D0*dy-9.D0/2.D0*dx*dy+3.D0/4.D0*(15.D0/2.D0*dx**2-3.D0/&
!                 2.D0)*dy-15.D0/8.D0*dy**3
!    Dhelp(9,1) = -3.D0*dx
!    Dhelp(10,1) = 3*(-1.D0/2.D0+3.D0/2.D0*dy**2)*dx
!
!    Dhelp(1,2) = -3.D0/4.D0+3.D0/2.D0*dy+3.D0/4.D0*dx**2
!    Dhelp(2,2) = -3.D0/2.D0*dx*dy
!    Dhelp(3,2) = 3.D0/4.D0+3.D0/2.D0*dy-3.D0/4.D0*dx**2
!    Dhelp(4,2) = 3.D0/2.D0*dx*dy
!    Dhelp(5,2) = -15.D0/8.D0*dx+9.D0/2.D0*dx*dy+15.D0/8.D0*dx**3-3.D0/4.D0*dx*&
!                 (15.D0/2.D0*dy**2-3.D0/2.D0)
!    Dhelp(6,2) = -3.D0/8.D0*dx-3.D0/4.D0+9.D0/4.D0*dx**2+15.D0/8.D0*dx**3-3.D0/&
!                 4.D0*dx*(15.D0/2.D0*dy**2-3.D0/2.D0)
!    Dhelp(7,2) = -15.D0/8.D0*dx-9.D0/2.D0*dx*dy+15.D0/8.D0*dx**3-3.D0/4.D0*dx*&
!                 (15.D0/2.D0*dy**2-3.D0/2.D0)
!    Dhelp(8,2) = -3.D0/8.D0*dx+3.D0/4.D0-9.D0/4.D0*dx**2+15.D0/8.D0*dx**3-3.D0/&
!                 4.D0*dx*(15.D0/2.D0*dy**2-3.D0/2.D0)
!    Dhelp(9,2) = -3*dy
!    Dhelp(10,2) = 3*(-1.D0/2.D0+3.D0/2.D0*dx**2)*dy

    ! x-derivatives on current element
!    if (Bder(DER_DERIV_X)) then
      Dbas( 1,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(1,1) - Djac(2)*Dhelp(1,2))
      Dbas( 2,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(2,1) - Djac(2)*Dhelp(2,2))
      Dbas( 3,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(3,1) - Djac(2)*Dhelp(3,2))
      Dbas( 4,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(4,1) - Djac(2)*Dhelp(4,2))
      Dbas( 5,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(5,1) - Djac(2)*Dhelp(5,2))*d5
      Dbas( 6,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(6,1) - Djac(2)*Dhelp(6,2))*d6
      Dbas( 7,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(7,1) - Djac(2)*Dhelp(7,2))*d7
      Dbas( 8,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(8,1) - Djac(2)*Dhelp(8,2))*d8
      Dbas( 9,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(9,1) - Djac(2)*Dhelp(9,2))
      Dbas(10,DER_DERIV_X) =  dxj*(Djac(4)*Dhelp(10,1) - Djac(2)*Dhelp(10,2))
!    endif

    ! y-derivatives on current element
!    if (Bder(DER_DERIV_Y)) then
      Dbas( 1,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(1,1) - Djac(1)*Dhelp(1,2))
      Dbas( 2,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(2,1) - Djac(1)*Dhelp(2,2))
      Dbas( 3,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(3,1) - Djac(1)*Dhelp(3,2))
      Dbas( 4,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(4,1) - Djac(1)*Dhelp(4,2))
      Dbas( 5,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(5,1) - Djac(1)*Dhelp(5,2))*d5
      Dbas( 6,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(6,1) - Djac(1)*Dhelp(6,2))*d6
      Dbas( 7,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(7,1) - Djac(1)*Dhelp(7,2))*d7
      Dbas( 8,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(8,1) - Djac(1)*Dhelp(8,2))*d8
      Dbas( 9,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(9,1) - Djac(1)*Dhelp(9,2))
      Dbas(10,DER_DERIV_Y) = -dxj*(Djac(3)*Dhelp(10,1) - Djac(1)*Dhelp(10,2))
!    endif
!  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_EB50_mult (celement, Dcoords, itwistIndex, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2TB.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges.
  integer(I32), intent(in) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(10,2,npoints) :: Dhelp

  ! Scaling factors for basis functions 5..8
  real(DP) :: d5,d6,d7,d8

  ! Auxiliary variables
  real(DP), dimension(npoints) :: dxj
  real(DP) :: dx,dy

  integer :: i   ! point counter

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: Q8 = 6.25_DP  ! = 25/4
  real(DP), parameter :: Q9 = 37.5_DP  ! = 75/2
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Get the twist indices for this element
  ! itwistIndex is a bitfield. Each bit specifies the orientation of an edge.
  ! We use the bit to calculate a "1.0" if the edge has positive orientation
  ! and "-1.0" if it has negative orientation.
  d5 = real(1-iand(int(ishft(itwistIndex, 1)),2),DP)
  d6 = real(1-iand(int(ishft(itwistIndex, 0)),2),DP)
  d7 = real(1-iand(int(ishft(itwistIndex,-1)),2),DP)
  d8 = real(1-iand(int(ishft(itwistIndex,-2)),2),DP)

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dbas(1,DER_FUNC,i) =  Q3*dy*( dx**2 + dy - P1) - Q1
      Dbas(2,DER_FUNC,i) =  Q3*dx*(-dy**2 + dx + P1) - Q1
      Dbas(3,DER_FUNC,i) =  Q3*dy*(-dx**2 + dy + P1) - Q1
      Dbas(4,DER_FUNC,i) =  Q3*dx*( dy**2 + dx - P1) - Q1
      Dbas(5,DER_FUNC,i) = -Q2*dx*(dy*(dy*( P4*dy - P5) - P4*dx**2 + P2) + P2)*d5
      Dbas(6,DER_FUNC,i) = -Q2*dy*(dx*(dx*(-P4*dx - P5) + P4*dy**2 - P2) + P2)*d6
      Dbas(7,DER_FUNC,i) = -Q2*dx*(dy*(dy*( P4*dy + P5) - P4*dx**2 + P2) - P2)*d7
      Dbas(8,DER_FUNC,i) = -Q2*dy*(dx*(dx*(-P4*dx + P5) + P4*dy**2 - P2) - P2)*d8
      Dbas(9,DER_FUNC,i) = -Q4*(dx**2 + dy**2) + P2
      Dbas(10,DER_FUNC,i)=  Q8*(P3*dx**2 - P1)*(P3*dy**2 - P1)
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN

    !x- and y-derivatives on reference element
    do i=1,npoints
      dxj(i) = 1.0_DP / Ddetj(i)
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dhelp(1,1,i) =  Q4*dx*dy
      Dhelp(2,1,i) = -Q3*dy**2 + Q4*dx + Q3
      Dhelp(3,1,i) = -Q4*dx*dy
      Dhelp(4,1,i) =  Q3*dy**2 + Q4*dx - Q3
      Dhelp(5,1,i) =  dy*(dy*(-Q6*dy + Q5) + Q7*dx**2 - Q3) - Q3
      Dhelp(6,1,i) = -Q2*dy*(dx*(-P7*dx - P6) + P4*dy**2 - P2)
      Dhelp(7,1,i) =  dy*(dy*(-Q6*dy - Q5) + Q7*dx**2 - Q3) + Q3
      Dhelp(8,1,i) = -Q2*dy*(dx*(-P7*dx + P6) + P4*dy**2 - P2)
      Dhelp(9,1,i) = -P3*dx
      Dhelp(10,1,i)=  Q9*dx*(P3*dy**2 - P1)
      Dhelp(1,2,i) =  Q3*dx**2 + Q4*dy - Q3
      Dhelp(2,2,i) = -Q4*dx*dy
      Dhelp(3,2,i) = -Q3*dx**2 + Q4*dy + Q3
      Dhelp(4,2,i) =  Q4*dx*dy
      Dhelp(5,2,i) = -Q2*dx*(dy*( P7*dy - P6) - P4*dx**2 + P2)
      Dhelp(6,2,i) =  dx*(dx*( Q6*dx + Q5) - Q7*dy**2 + Q3) - Q3
      Dhelp(7,2,i) =  Q2*dx*(dy*(-P7*dy - P6) + P4*dx**2 - P2)
      Dhelp(8,2,i) =  dx*(dx*( Q6*dx - Q5) - Q7*dy**2 + Q3) + Q3
      Dhelp(9,2,i) = -P3*dy
      Dhelp(10,2,i)=  Q9*dy*(P3*dx**2 - P1)
    end do

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        Dbas(1,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(1,1,i) &
                                         - Djac(2,i)*Dhelp(1,2,i))
        Dbas(2,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(2,1,i) &
                                         - Djac(2,i)*Dhelp(2,2,i))
        Dbas(3,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(3,1,i) &
                                         - Djac(2,i)*Dhelp(3,2,i))
        Dbas(4,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(4,1,i) &
                                         - Djac(2,i)*Dhelp(4,2,i))
        Dbas(5,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(5,1,i) &
                                         - Djac(2,i)*Dhelp(5,2,i))*d5
        Dbas(6,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(6,1,i) &
                                         - Djac(2,i)*Dhelp(6,2,i))*d6
        Dbas(7,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(7,1,i) &
                                         - Djac(2,i)*Dhelp(7,2,i))*d7
        Dbas(8,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(8,1,i) &
                                         - Djac(2,i)*Dhelp(8,2,i))*d8
        Dbas(9,DER_DERIV_X,i) =  dxj(i)*(Djac(4,i)*Dhelp(9,1,i) &
                                         - Djac(2,i)*Dhelp(9,2,i))
        Dbas(10,DER_DERIV_X,i)=  dxj(i)*(Djac(4,i)*Dhelp(10,1,i) &
                                         - Djac(2,i)*Dhelp(10,2,i))
!      END DO
!    ENDIF

    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(1,1,i) &
                                         - Djac(1,i)*Dhelp(1,2,i))
        Dbas(2,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(2,1,i) &
                                         - Djac(1,i)*Dhelp(2,2,i))
        Dbas(3,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(3,1,i) &
                                         - Djac(1,i)*Dhelp(3,2,i))
        Dbas(4,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(4,1,i) &
                                         - Djac(1,i)*Dhelp(4,2,i))
        Dbas(5,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(5,1,i) &
                                         - Djac(1,i)*Dhelp(5,2,i))*d5
        Dbas(6,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(6,1,i) &
                                         - Djac(1,i)*Dhelp(6,2,i))*d6
        Dbas(7,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(7,1,i) &
                                         - Djac(1,i)*Dhelp(7,2,i))*d7
        Dbas(8,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(8,1,i) &
                                         - Djac(1,i)*Dhelp(8,2,i))*d8
        Dbas(9,DER_DERIV_Y,i) = -dxj(i)*(Djac(3,i)*Dhelp(9,1,i) &
                                         - Djac(1,i)*Dhelp(9,2,i))
        Dbas(10,DER_DERIV_Y,i)= -dxj(i)*(Djac(3,i)*Dhelp(10,1,i) &
                                         - Djac(1,i)*Dhelp(10,2,i))
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_EB50_sim (celement, Dcoords, ItwistIndex, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2TB.
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

  ! List of twist indices. For every edge/face on every cell, the twist
  ! index defines the orientation of the edge/face.
  ! Array with DIMENSION(1:NVE/NVA,nelements)
  integer(I32), dimension(:), intent(in) :: ItwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(10,2,npoints) :: Dhelp

  ! Scaling factors for basis functions 5..8
  real(DP) :: d5,d6,d7,d8

  ! Auxiliary variables
  real(DP), dimension(npoints) :: dxj
  real(DP) :: dx,dy

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: Q8 = 6.25_DP  ! = 25/4
  real(DP), parameter :: Q9 = 37.5_DP  ! = 75/2
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i,dx,dy,d5,d6,d7,d8) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      ! Get the twist indices for this element.
      ! ItwistIndex(.) is a bitfield. Each bit specifies the orientation of an edge.
      ! We use the bit to calculate a "1.0" if the edge has positive orientation
      ! and "-1.0" if it has negative orientation.
      d5 = real(1-iand(int(ishft(ItwistIndex(j), 1)),2),DP)
      d6 = real(1-iand(int(ishft(ItwistIndex(j), 0)),2),DP)
      d7 = real(1-iand(int(ishft(ItwistIndex(j),-1)),2),DP)
      d8 = real(1-iand(int(ishft(ItwistIndex(j),-2)),2),DP)

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dbas(1,DER_FUNC,i,j) =  Q3*dy*( dx**2 + dy - P1) - Q1
        Dbas(2,DER_FUNC,i,j) =  Q3*dx*(-dy**2 + dx + P1) - Q1
        Dbas(3,DER_FUNC,i,j) =  Q3*dy*(-dx**2 + dy + P1) - Q1
        Dbas(4,DER_FUNC,i,j) =  Q3*dx*( dy**2 + dx - P1) - Q1
        Dbas(5,DER_FUNC,i,j) = -Q2*dx*(dy*(dy*( P4*dy - P5) - P4*dx**2 + P2) + P2)*d5
        Dbas(6,DER_FUNC,i,j) = -Q2*dy*(dx*(dx*(-P4*dx - P5) + P4*dy**2 - P2) + P2)*d6
        Dbas(7,DER_FUNC,i,j) = -Q2*dx*(dy*(dy*( P4*dy + P5) - P4*dx**2 + P2) - P2)*d7
        Dbas(8,DER_FUNC,i,j) = -Q2*dy*(dx*(dx*(-P4*dx + P5) + P4*dy**2 - P2) - P2)*d8
        Dbas(9,DER_FUNC,i,j) = -Q4*(dx**2 + dy**2) + P2
        Dbas(10,DER_FUNC,i,j)=  Q8*(P3*dx**2 - P1)*(P3*dy**2 - P1)
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dx,dy,d5,d6,d7,d8,dxj,Dhelp) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      ! Get the twist indices for this element
      ! ItwistIndex(.) is a bitfield. Each bit specifies the orientation of an edge.
      ! We use the bit to calculate a "1.0" if the edge has positive orientation
      ! and "-1.0" if it has negative orientation.
      d5 = real(1-iand(int(ishft(ItwistIndex(j), 1)),2),DP)
      d6 = real(1-iand(int(ishft(ItwistIndex(j), 0)),2),DP)
      d7 = real(1-iand(int(ishft(ItwistIndex(j),-1)),2),DP)
      d8 = real(1-iand(int(ishft(ItwistIndex(j),-2)),2),DP)

      !x- and y-derivatives on reference element
      do i=1,npoints
        dxj(i) = 1.0_DP / Ddetj(i,j)
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dhelp(1,1,i) =  Q4*dx*dy
        Dhelp(2,1,i) = -Q3*dy**2 + Q4*dx + Q3
        Dhelp(3,1,i) = -Q4*dx*dy
        Dhelp(4,1,i) =  Q3*dy**2 + Q4*dx - Q3
        Dhelp(5,1,i) =  dy*(dy*(-Q6*dy + Q5) + Q7*dx**2 - Q3) - Q3
        Dhelp(6,1,i) = -Q2*dy*(dx*(-P7*dx - P6) + P4*dy**2 - P2)
        Dhelp(7,1,i) =  dy*(dy*(-Q6*dy - Q5) + Q7*dx**2 - Q3) + Q3
        Dhelp(8,1,i) = -Q2*dy*(dx*(-P7*dx + P6) + P4*dy**2 - P2)
        Dhelp(9,1,i) = -P3*dx
        Dhelp(10,1,i)=  Q9*dx*(P3*dy**2 - P1)
        Dhelp(1,2,i) =  Q3*dx**2 + Q4*dy - Q3
        Dhelp(2,2,i) = -Q4*dx*dy
        Dhelp(3,2,i) = -Q3*dx**2 + Q4*dy + Q3
        Dhelp(4,2,i) =  Q4*dx*dy
        Dhelp(5,2,i) = -Q2*dx*(dy*( P7*dy - P6) - P4*dx**2 + P2)
        Dhelp(6,2,i) =  dx*(dx*( Q6*dx + Q5) - Q7*dy**2 + Q3) - Q3
        Dhelp(7,2,i) =  Q2*dx*(dy*(-P7*dy - P6) + P4*dx**2 - P2)
        Dhelp(8,2,i) =  dx*(dx*( Q6*dx - Q5) - Q7*dy**2 + Q3) + Q3
        Dhelp(9,2,i) = -P3*dy
        Dhelp(10,2,i)=  Q9*dy*(P3*dx**2 - P1)
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(1,1,i) &
                                           - Djac(2,i,j)*Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(2,1,i) &
                                           - Djac(2,i,j)*Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(3,1,i) &
                                           - Djac(2,i,j)*Dhelp(3,2,i))
          Dbas(4,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(4,1,i) &
                                           - Djac(2,i,j)*Dhelp(4,2,i))
          Dbas(5,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(5,1,i) &
                                           - Djac(2,i,j)*Dhelp(5,2,i))*d5
          Dbas(6,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(6,1,i) &
                                           - Djac(2,i,j)*Dhelp(6,2,i))*d6
          Dbas(7,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(7,1,i) &
                                           - Djac(2,i,j)*Dhelp(7,2,i))*d7
          Dbas(8,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(8,1,i) &
                                           - Djac(2,i,j)*Dhelp(8,2,i))*d8
          Dbas(9,DER_DERIV_X,i,j) =  dxj(i)*(Djac(4,i,j)*Dhelp(9,1,i) &
                                           - Djac(2,i,j)*Dhelp(9,2,i))
          Dbas(10,DER_DERIV_X,i,j)=  dxj(i)*(Djac(4,i,j)*Dhelp(10,1,i) &
                                           - Djac(2,i,j)*Dhelp(10,2,i))
        !end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        !do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(1,1,i) &
                                           - Djac(1,i,j)*Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(2,1,i) &
                                           - Djac(1,i,j)*Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(3,1,i) &
                                           - Djac(1,i,j)*Dhelp(3,2,i))
          Dbas(4,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(4,1,i) &
                                           - Djac(1,i,j)*Dhelp(4,2,i))
          Dbas(5,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(5,1,i) &
                                           - Djac(1,i,j)*Dhelp(5,2,i))*d5
          Dbas(6,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(6,1,i) &
                                           - Djac(1,i,j)*Dhelp(6,2,i))*d6
          Dbas(7,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(7,1,i) &
                                           - Djac(1,i,j)*Dhelp(7,2,i))*d7
          Dbas(8,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(8,1,i) &
                                           - Djac(1,i,j)*Dhelp(8,2,i))*d8
          Dbas(9,DER_DERIV_Y,i,j) = -dxj(i)*(Djac(3,i,j)*Dhelp(9,1,i) &
                                           - Djac(1,i,j)*Dhelp(9,2,i))
          Dbas(10,DER_DERIV_Y,i,j)= -dxj(i)*(Djac(3,i,j)*Dhelp(10,1,i) &
                                           - Djac(1,i,j)*Dhelp(10,2,i))
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric DG_T0_2D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_DG_T0_2D (celement, Dcoords, Djac, ddetj, Bder, &
                                 Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T0_2D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  ! Clear the output array
  !Dbas = 0.0_DP

  ! DG_T0_2D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T0_2D_mult (celement, Dcoords, Djac, Ddetj, &
                                      Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T0_2D.
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
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints

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

  ! Clear the output array
  !Dbas = 0.0_DP

  ! DG_T0_2D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T0_2D_sim (celement, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T0_2D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
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
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  !</input>

  !<output>

  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th
  !   basis function of the finite element in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas

  !</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP

  ! DG_T0_2D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  end subroutine

!**************************************************************************
! Element subroutines for parametric DG_T1_2D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_DG_T1_2D (celement, Dcoords, Djac, ddetj, Bder, &
                                 Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T1_2D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D) :: Dhelp
  real(DP) :: dx,dy

  real(DP) :: dxj !auxiliary variable

  ! The DG_T1_2D element is specified by three polynomials on the reference element.
  ! These three polynomials are:
  !
  !  P1(X,Y) = 1
  !  P2(X,Y) = x
  !  P3(X,Y) = y
  !
  ! The first one is the constant function, which can be used to set the value of
  ! the solution in the centroid of the element and the second and third have
  ! vanishing mean-values on the reference element and can be used to set the
  ! x and y derivatives in the centroid of the element while not changing the
  ! mean-value in the element.
  ! So slope limiting is possible without affecting the conservation laws.

  ! Clear the output array
  !Dbas = 0.0_DP

  dx = Dpoint(1)
  dy = Dpoint(2)

  ! Remark: The DG_T1_2D-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 1.0_DP
    Dbas(2,DER_FUNC) = dx
    Dbas(3,DER_FUNC) = dy
!  endif

  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 1.0E0_DP / ddetj          !!!!!!!!!!!! 0.25E0_DP / ddetj

    ! x- and y-derivatives on reference element
    Dhelp(1,1) = 0.0_DP
    Dhelp(2,1) = 1.0_DP
    Dhelp(3,1) = 0.0_DP
    Dhelp(1,2) = 0.0_DP
    Dhelp(2,2) = 0.0_DP
    Dhelp(3,2) = 1.0_DP

    ! x-derivatives on current element
!    if (Bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(1,1) - Djac(2) * Dhelp(1,2))
      Dbas(2,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(2,1) - Djac(2) * Dhelp(2,2))
      Dbas(3,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(3,1) - Djac(2) * Dhelp(3,2))
!    endif

    ! y-derivatives on current element
!    if (Bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(1,1) - Djac(1) * Dhelp(1,2))
      Dbas(2,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(2,1) - Djac(1) * Dhelp(2,2))
      Dbas(3,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(3,1) - Djac(1) * Dhelp(3,2))
!    endif
!  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T1_2D_mult (celement, Dcoords, Djac, Ddetj, &
                                      Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T1_2D.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(3,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T1_2D-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC,i) = 1.0_DP
      Dbas(2,DER_FUNC,i) = Dpoints(1,i)
      Dbas(3,DER_FUNC,i) = Dpoints(2,i)
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj(:) = 1.0E0_DP / Ddetj(1:npoints)        !!!!!!!!!!Dxj(:) = 0.25E0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
      Dhelp(1,1,i) = 0.0_DP
      Dhelp(2,1,i) = 1.0_DP
      Dhelp(3,1,i) = 0.0_DP
      Dhelp(1,2,i) = 0.0_DP
      Dhelp(2,2,i) = 0.0_DP
      Dhelp(3,2,i) = 1.0_DP
    end do

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        Dbas(1,DER_DERIV_X,i) = dxj(i) * (Djac(4,i) * Dhelp(1,1,i) - Djac(2,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_X,i) = dxj(i) * (Djac(4,i) * Dhelp(2,1,i) - Djac(2,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_X,i) = dxj(i) * (Djac(4,i) * Dhelp(3,1,i) - Djac(2,i) * Dhelp(3,2,i))
!      END DO
!    ENDIF

    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = -dxj(i) * (Djac(3,i) * Dhelp(1,1,i) - Djac(1,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_Y,i) = -dxj(i) * (Djac(3,i) * Dhelp(2,1,i) - Djac(1,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_Y,i) = -dxj(i) * (Djac(3,i) * Dhelp(3,1,i) - Djac(1,i) * Dhelp(3,2,i))
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_DG_T1_2D_sim (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, nelements, &
                                Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T1_2D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(3,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T1_2D-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 1.0_DP
        Dbas(2,DER_FUNC,i,j) = Dpoints(1,i,j)
        Dbas(3,DER_FUNC,i,j) = Dpoints(2,i,j)
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj,Dhelp) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1.0E0_DP / Ddetj(1:npoints,j) !!!!!!!!! Dxj(:) = 0.25E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
        Dhelp(1,1,i) = 0.0_DP
        Dhelp(2,1,i) = 1.0_DP
        Dhelp(3,1,i) = 0.0_DP
        Dhelp(1,2,i) = 0.0_DP
        Dhelp(2,2,i) = 0.0_DP
        Dhelp(3,2,i) = 1.0_DP
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(1,1,i) &
                                    - Djac(2,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(2,1,i) &
                                    - Djac(2,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(3,1,i) &
                                    - Djac(2,i,j) * Dhelp(3,2,i))
        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(1,1,i) &
                                    - Djac(1,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(2,1,i) &
                                    - Djac(1,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(3,1,i) &
                                    - Djac(1,i,j) * Dhelp(3,2,i))
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine


!**************************************************************************
! Element subroutines for parametric DG_T2_2D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_DG_T2_2D (celement, Dcoords, Djac, ddetj, Bder, &
                                 Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T2_2D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,5) :: Dhelp
  real(DP) :: dx,dy

  integer :: idof
  real(DP) :: dxj,dxjs,dc1,dc2,dc3 !auxiliary variable

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q3 = 1.0_DP/3.0_DP
  real(DP), parameter :: Q4 = 0.25_DP
  real(DP), parameter :: Q6 = 1.0_DP/6.0_DP

  ! The DG_T2_2D element is specified by six polynomials on the reference element.
  ! These six polynomials are:
  !
  !  P1(X,Y) = 1
  !  P2(X,Y) = x
  !  P3(X,Y) = y
  !  P4(X,Y) = 1/2 (x^2 - 1/3)
  !  P5(X,Y) = 1/2 (y^2 - 1/3)
  !  P6(X,Y) = xy
  !
  ! The first one is the constant function, which can be used to set the value of
  ! the solution in the centroid of the element and the second and third have
  ! vanishing mean-values on the reference element and can be used to set the
  ! x and y derivatives in the centroid of the element while not changing the
  ! mean-value in the element.
  ! Analogous the fourth to sixth polynomial can be used to set the second
  ! derivative in xx, xy and yy direction in the centroid without affecting the
  ! first derivatives or mean-values.
  ! So slope limiting is possible without affecting the conservation laws.

  ! Clear the output array
  !Dbas = 0.0_DP

  dx = Dpoint(1)
  dy = Dpoint(2)

  !if function values are desired
  if (Bder(DER_FUNC)) then
    Dbas(1,DER_FUNC)= 1.0_DP
    Dbas(2,DER_FUNC)= dx
    Dbas(3,DER_FUNC)= dy
    Dbas(4,DER_FUNC)= Q2*dx*dx-Q6
    Dbas(5,DER_FUNC)= Q2*dy*dy-Q6
    Dbas(6,DER_FUNC)= dx*dy
  endif

  ! if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 1.0E0_DP / ddetj

    !x- and y-derivatives on reference element
    Dhelp(1,1)= 0.0_DP
    Dhelp(2,1)= 1.0_DP
    Dhelp(3,1)= 0.0_DP
    Dhelp(4,1)= dx
    Dhelp(5,1)= 0.0_DP
    Dhelp(6,1)= dy

    Dhelp(1,2)= 0.0_DP
    Dhelp(2,2)= 0.0_DP
    Dhelp(3,2)= 1.0_DP
    Dhelp(4,2)= 0.0_DP
    Dhelp(5,2)= dy
    Dhelp(6,2)= dx

    ! x-derivatives on current element
    if (Bder(DER_DERIV_X)) then
      do idof = 1,6
        Dbas(idof,DER_DERIV_X) = &
            dxj * (Djac(4) * Dhelp(idof,1) - Djac(2) * Dhelp(idof,2))
      end do
    endif

    ! y-derivatives on current element
    if (Bder(DER_DERIV_Y)) then
      do idof = 1,6
        Dbas(idof,DER_DERIV_Y) = &
            -dxj * (Djac(3) * Dhelp(idof,1) - Djac(1) * Dhelp(idof,2))
      end do
    endif
  endif

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then
    dxj = 1.0_DP / ddetj

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, our finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]
    !

    dxjs = dxj*dxj

    !xx-, xy and yy-derivatives on reference element
    Dhelp(1,3) = 0.0_DP
    Dhelp(2,3) = 0.0_DP
    Dhelp(3,3) = 0.0_DP
    Dhelp(4,3) = 1.0_DP
    Dhelp(5,3) = 0.0_DP
    Dhelp(6,3) = 0.0_DP

    Dhelp(1,4) = 0.0_DP
    Dhelp(2,4) = 0.0_DP
    Dhelp(3,4) = 0.0_DP
    Dhelp(4,4) = 0.0_DP
    Dhelp(5,4) = 0.0_DP
    Dhelp(6,4) = 1.0_DP

    Dhelp(1,5) = 0.0_DP
    Dhelp(2,5) = 0.0_DP
    Dhelp(3,5) = 0.0_DP
    Dhelp(4,5) = 0.0_DP
    Dhelp(5,5) = 1.0_DP
    Dhelp(6,5) = 0.0_DP

    ! WARNING: NOT TESTED!!!

    ! xx-derivatives on current element
    if (Bder(DER_DERIV2D_XX)) then
      dc1 = dxjs * Djac(4)**2
      dc2 = dxjs * ( -2.0_DP * Djac(4) * Djac(2) )
      dc3 = dxjs * Djac(2)**2
      do idof = 1,6
        Dbas(idof,DER_DERIV2D_XX) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif

    ! xy-derivatives on current element
    if (Bder(DER_DERIV2D_XY)) then
      dc1 = - dxjs * Djac(4) * Djac(3)
      dc2 = dxjs * ( Djac(4) * Djac(1) + Djac(3) * Djac(2) )
      dc3 = - dxjs * Djac(1) * Djac(4)
      do idof = 1,6
        Dbas(idof,DER_DERIV2D_XY) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif

    ! yy-derivatives on current element
    if (Bder(DER_DERIV2D_YY)) then
      dc1 = dxjs * Djac(3)**2
      dc2 = dxjs * ( -2.0_DP * Djac(3) * Djac(1) )
      dc3 = dxjs * Djac(1)**2
      do idof = 1,6
        Dbas(idof,DER_DERIV2D_YY) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif
  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T2_2D_mult (celement, Dcoords, Djac, Ddetj, &
                                      Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T2_2D.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,5,npoints) :: Dhelp
  real(DP) :: dx,dy
  real(DP),dimension(npoints) :: Dxj,Dxjs,Dc1,Dc2,Dc3 !auxiliary variable

  integer :: i,idof   ! point counter

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q3 = 1.0_DP/3.0_DP
  real(DP), parameter :: Q4 = 0.25_DP
  real(DP), parameter :: Q6 = 1.0_DP/6.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T1_2D-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dbas(1,DER_FUNC,i)= 1.0_DP
      Dbas(2,DER_FUNC,i)= dx
      Dbas(3,DER_FUNC,i)= dy
      Dbas(4,DER_FUNC,i)= Q2*dx*dx-Q6
      Dbas(5,DER_FUNC,i)= Q2*dy*dy-Q6
      Dbas(6,DER_FUNC,i)= dx*dy
    end do
  !ENDIF

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    Dxj(:) = 1.0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)

      !x- and y-derivatives on reference element
      Dhelp(1,1,i)= 0.0_DP
      Dhelp(2,1,i)= 1.0_DP
      Dhelp(3,1,i)= 0.0_DP
      Dhelp(4,1,i)= dx
      Dhelp(5,1,i)= 0.0_DP
      Dhelp(6,1,i)= dy

      Dhelp(1,2,i)= 0.0_DP
      Dhelp(2,2,i)= 0.0_DP
      Dhelp(3,2,i)= 1.0_DP
      Dhelp(4,2,i)= 0.0_DP
      Dhelp(5,2,i)= dy
      Dhelp(6,2,i)= dx
    end do

    ! x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        do idof = 1,6
          Dbas(idof,DER_DERIV_X,i) = &
              Dxj(i) * (Djac(4,i) * Dhelp(idof,1,i) - Djac(2,i) * Dhelp(idof,2,i))
        end do
!      END DO
!    ENDIF

    ! y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        do idof = 1,6
          Dbas(idof,DER_DERIV_Y,i) = &
              -Dxj(i) * (Djac(3,i) * Dhelp(idof,1,i) - Djac(1,i) * Dhelp(idof,2,i))
        end do
      end do
!    ENDIF
  endif

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]

    !x- and y-derivatives on reference element
    Dxj(:) = 1.0E0_DP / Ddetj(1:npoints)
    Dxjs = Dxj*Dxj

    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)

      !xx-, xy and yy-derivatives on reference element
          Dhelp(1,3,i) = 0.0_DP
          Dhelp(2,3,i) = 0.0_DP
          Dhelp(3,3,i) = 0.0_DP
          Dhelp(4,3,i) = 1.0_DP
          Dhelp(5,3,i) = 0.0_DP
          Dhelp(6,3,i) = 0.0_DP

          Dhelp(1,4,i) = 0.0_DP
          Dhelp(2,4,i) = 0.0_DP
          Dhelp(3,4,i) = 0.0_DP
          Dhelp(4,4,i) = 0.0_DP
          Dhelp(5,4,i) = 0.0_DP
          Dhelp(6,4,i) = 1.0_DP

          Dhelp(1,5,i) = 0.0_DP
          Dhelp(2,5,i) = 0.0_DP
          Dhelp(3,5,i) = 0.0_DP
          Dhelp(4,5,i) = 0.0_DP
          Dhelp(5,5,i) = 1.0_DP
          Dhelp(6,5,i) = 0.0_DP

    end do

    ! WARNING: NOT TESTED!!!

    ! xx-derivatives on current element
    !if (Bder(DER_DERIV_XX)) then
      do i=1,npoints
        Dc1(i) = Dxjs(i) * Djac(4,i)**2
        Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(4,i) * Djac(2,i) )
        Dc3(i) = Dxjs(i) * Djac(2,i)**2
      end do

      do i=1,npoints
        do idof = 1,6
          Dbas(idof,DER_DERIV_XX,i) = &
              Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
        end do
      end do
    !endif

    ! xy-derivatives on current element
    !if (Bder(DER_DERIV_XY)) then
      do i=1,npoints
        Dc1(i) = - Dxjs(i) * Djac(4,i) * Djac(3,i)
        Dc2(i) = Dxjs(i) * ( Djac(4,i) * Djac(1,i) + Djac(3,i) * Djac(2,i) )
        Dc3(i) = - Dxjs(i) * Djac(1,i) * Djac(4,i)
      end do

      do i=1,npoints
        do idof = 1,6
          Dbas(idof,DER_DERIV_XY,i) = &
              Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
        end do
      end do
    !endif

    ! yy-derivatives on current element
    !if (Bder(DER_DERIV_YY)) then
      do i=1,npoints
        Dc1(i) = Dxjs(i) * Djac(3,i)**2
        Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(3,i) * Djac(1,i) )
        Dc3(i) = Dxjs(i) * Djac(1,i)**2
      end do

      do idof = 1,6
        Dbas(idof,DER_DERIV_YY,i) = &
            Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
      end do
    !endif
  endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_DG_T2_2D_sim (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, nelements, &
                                Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T2_2D.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,5,npoints) :: Dhelp
  real(DP) :: dx,dy
  real(DP),dimension(npoints) :: Dxj,Dxjs,Dc1,Dc2,Dc3 !auxiliary variables

  integer :: i   ! point counter
  integer :: j   ! element counter
  integer :: idof

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q3 = 1.0_DP/3.0_DP
  real(DP), parameter :: Q4 = 0.25_DP
  real(DP), parameter :: Q6 = 1.0_DP/6.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T1_2D-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i,dx,dy) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dbas(1,DER_FUNC,i,j)= 1.0_DP
        Dbas(2,DER_FUNC,i,j)= dx
        Dbas(3,DER_FUNC,i,j)= dy
        Dbas(4,DER_FUNC,i,j)= Q2*dx*dx-Q6
        Dbas(5,DER_FUNC,i,j)= Q2*dy*dy-Q6
        Dbas(6,DER_FUNC,i,j)= dx*dy
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,Dxj,dx,dy,Dhelp,idof) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1.0E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)

        !x- and y-derivatives on reference element
        Dhelp(1,1,i)= 0.0_DP
        Dhelp(2,1,i)= 1.0_DP
        Dhelp(3,1,i)= 0.0_DP
        Dhelp(4,1,i)= dx
        Dhelp(5,1,i)= 0.0_DP
        Dhelp(6,1,i)= dy

        Dhelp(1,2,i)= 0.0_DP
        Dhelp(2,2,i)= 0.0_DP
        Dhelp(3,2,i)= 1.0_DP
        Dhelp(4,2,i)= 0.0_DP
        Dhelp(5,2,i)= dy
        Dhelp(6,2,i)= dx
      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          do idof = 1,6
            Dbas(idof,DER_DERIV_X,i,j) = &
                Dxj(i) * (Djac(4,i,j) * Dhelp(idof,1,i) &
                          - Djac(2,i,j) * Dhelp(idof,2,i))
!          end do
!        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
!        do i=1,npoints
!          do idof = 1,6
            Dbas(idof,DER_DERIV_Y,i,j) = &
                -Dxj(i) * (Djac(3,i,j) * Dhelp(idof,1,i) &
                - Djac(1,i,j) * Dhelp(idof,2,i))
          end do
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]

    !x- and y-derivatives on reference element
    !$omp parallel do default(shared)&
    !$omp private(Dc1,Dc2,Dc3,Dhelp,Dxj,Dxjs,dx,dy,i,idof)&
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1.0E0_DP / Ddetj(1:npoints,j)
      Dxjs = Dxj*Dxj

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)

        !xx-, xy and yy-derivatives on reference element
        Dhelp(1,3,i) = 0.0_DP
        Dhelp(2,3,i) = 0.0_DP
        Dhelp(3,3,i) = 0.0_DP
        Dhelp(4,3,i) = 1.0_DP
        Dhelp(5,3,i) = 0.0_DP
        Dhelp(6,3,i) = 0.0_DP

        Dhelp(1,4,i) = 0.0_DP
        Dhelp(2,4,i) = 0.0_DP
        Dhelp(3,4,i) = 0.0_DP
        Dhelp(4,4,i) = 0.0_DP
        Dhelp(5,4,i) = 0.0_DP
        Dhelp(6,4,i) = 1.0_DP

        Dhelp(1,5,i) = 0.0_DP
        Dhelp(2,5,i) = 0.0_DP
        Dhelp(3,5,i) = 0.0_DP
        Dhelp(4,5,i) = 0.0_DP
        Dhelp(5,5,i) = 1.0_DP
        Dhelp(6,5,i) = 0.0_DP

      end do

      ! WARNING: NOT TESTED!!!

      ! xx-derivatives on current element
      !if (Bder(DER_DERIV_XX)) then
        do i=1,npoints
          Dc1(i) = Dxjs(i) * Djac(4,i,j)**2
          Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(4,i,j) * Djac(2,i,j) )
          Dc3(i) = Dxjs(i) * Djac(2,i,j)**2
        end do

        do i=1,npoints
          do idof = 1,6
            Dbas(idof,DER_DERIV_XX,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

      ! xy-derivatives on current element
      !if (Bder(DER_DERIV_XY)) then
        do i=1,npoints
          Dc1(i) = - Dxjs(i) * Djac(4,i,j) * Djac(3,i,j)
          Dc2(i) = Dxjs(i) * ( Djac(4,i,j) * Djac(1,i,j) + Djac(3,i,j) * Djac(2,i,j) )
          Dc3(i) = - Dxjs(i) * Djac(1,i,j) * Djac(4,i,j)
        end do

        do i=1,npoints
          do idof = 1,6
            Dbas(idof,DER_DERIV_XY,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

      ! yy-derivatives on current element
      !if (Bder(DER_DERIV_YY)) then
        do i=1,npoints
          Dc1(i) = Dxjs(i) * Djac(3,i,j)**2
          Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(3,i,j) * Djac(1,i,j) )
          Dc3(i) = Dxjs(i) * Djac(1,i,j)**2
        end do

        do i=1,npoints
          do idof = 1,6
            Dbas(idof,DER_DERIV_YY,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

    end do
    !$omp end parallel do

  end if

  end subroutine


!**************************************************************************
! Element subroutines for parametric DG_T3_2D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_DG_T3_2D (celement, Dcoords, Djac, ddetj, Bder, &
                                 Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T3_2D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(in) :: Dpoint
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

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(10,5) :: Dhelp
  real(DP) :: dx,dy

  integer :: idof
  real(DP) :: dxj,dxjs,dc1,dc2,dc3 !auxiliary variable

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q3 = 1.0_DP/3.0_DP
  real(DP), parameter :: Q4 = 0.25_DP
  real(DP), parameter :: Q6 = 1.0_DP/6.0_DP
  real(DP), parameter :: Q48 = 1.0_DP/48.0_DP

  ! The DG_T2_2D element is specified by six polynomials on the reference element.
  ! These six polynomials are:
  !
  !  P1(X,Y) = 1
  !  P2(X,Y) = x
  !  P3(X,Y) = y
  !  P4(X,Y) = 1/2 (x^2 - 1/3)
  !  P5(X,Y) = 1/2 (y^2 - 1/3)
  !  P6(X,Y) = xy
  !  P7(X,Y) = 1/48 x^3
  !  P8(X,Y) = 1/48 x^2 y
  !  P9(X,Y) = 1/48 x y^2
  !  P10(X,Y) = 1/48 y^3
  !
  ! The first one is the constant function, which can be used to set the value of
  ! the solution in the centroid of the element and the second and third have
  ! vanishing mean-values on the reference element and can be used to set the
  ! x and y derivatives in the centroid of the element while not changing the
  ! mean-value in the element.
  ! Analogous the fourth to sixth polynomial can be used to set the second
  ! derivative in xx, xy and yy direction in the centroid without affecting the
  ! first derivatives or mean-values.
  ! So slope limiting is possible without affecting the conservation laws.

  ! Clear the output array
  !Dbas = 0.0_DP

  dx = Dpoint(1)
  dy = Dpoint(2)

  !if function values are desired
  if (Bder(DER_FUNC)) then
    Dbas(1,DER_FUNC)= 1.0_DP
    Dbas(2,DER_FUNC)= dx
    Dbas(3,DER_FUNC)= dy
    Dbas(4,DER_FUNC)= Q2*dx*dx-Q6
    Dbas(5,DER_FUNC)= Q2*dy*dy-Q6
    Dbas(6,DER_FUNC)= dx*dy
    Dbas(7,DER_FUNC)= Q48*dx*dx*dx
    Dbas(8,DER_FUNC)= Q48*dx*dx*dy
    Dbas(9,DER_FUNC)= Q48*dx*dy*dy
    Dbas(10,DER_FUNC)= Q48*dy*dy*dy
  endif

  ! if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 1.0E0_DP / ddetj

    !x- and y-derivatives on reference element
    Dhelp(1,1)= 0.0_DP
    Dhelp(2,1)= 1.0_DP
    Dhelp(3,1)= 0.0_DP
    Dhelp(4,1)= dx
    Dhelp(5,1)= 0.0_DP
    Dhelp(6,1)= dy
    Dhelp(7,1)= 1.0_dp/16.0_dp*dx*dx
    Dhelp(8,1)= 1.0_dp/24.0_dp*dx*dy
    Dhelp(9,1)= 1.0_dp/48.0_dp*dy*dy
    Dhelp(10,1)= 0.0_dp

    Dhelp(1,2)= 0.0_DP
    Dhelp(2,2)= 0.0_DP
    Dhelp(3,2)= 1.0_DP
    Dhelp(4,2)= 0.0_DP
    Dhelp(5,2)= dy
    Dhelp(6,2)= dx
    Dhelp(7,2)= 0.0_dp
    Dhelp(8,2)= 1.0_dp/48.0_dp*dx*dx
    Dhelp(9,2)= 1.0_dp/24.0_dp*dx*dy
    Dhelp(10,2)= 1.0_dp/16.0_dp*dy*dy

    ! x-derivatives on current element
    if (Bder(DER_DERIV_X)) then
      do idof = 1,10
        Dbas(idof,DER_DERIV_X) = &
            dxj * (Djac(4) * Dhelp(idof,1) - Djac(2) * Dhelp(idof,2))
      end do
    endif

    ! y-derivatives on current element
    if (Bder(DER_DERIV_Y)) then
      do idof = 1,10
        Dbas(idof,DER_DERIV_Y) = &
            -dxj * (Djac(3) * Dhelp(idof,1) - Djac(1) * Dhelp(idof,2))
      end do
    endif
  endif

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then
    dxj = 1.0_DP / ddetj

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, our finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]
    !

    dxjs = dxj*dxj

    !xx-, xy and yy-derivatives on reference element
    Dhelp(1,3) = 0.0_DP
    Dhelp(2,3) = 0.0_DP
    Dhelp(3,3) = 0.0_DP
    Dhelp(4,3) = 1.0_DP
    Dhelp(5,3) = 0.0_DP
    Dhelp(6,3) = 0.0_DP
    Dhelp(7,3) = 1.0_DP/8.0_dp*dx
    Dhelp(8,3) = 1.0_DP/24.0_dp*dy
    Dhelp(9,3) = 0.0_DP
    Dhelp(10,3) = 0.0_DP

    Dhelp(1,4) = 0.0_DP
    Dhelp(2,4) = 0.0_DP
    Dhelp(3,4) = 0.0_DP
    Dhelp(4,4) = 0.0_DP
    Dhelp(5,4) = 0.0_DP
    Dhelp(6,4) = 1.0_DP
    Dhelp(7,4) = 0.0_dp
    Dhelp(8,4) = 1.0_DP/24.0_dp*dx
    Dhelp(9,4) = 1.0_DP/24.0_dp*dy
    Dhelp(10,4) = 0.0_DP

    Dhelp(1,5) = 0.0_DP
    Dhelp(2,5) = 0.0_DP
    Dhelp(3,5) = 0.0_DP
    Dhelp(4,5) = 0.0_DP
    Dhelp(5,5) = 1.0_DP
    Dhelp(6,5) = 0.0_DP
    Dhelp(7,5) = 0.0_dp
    Dhelp(8,5) = 0.0_dp
    Dhelp(9,5) = 1.0_DP/24.0_dp*dx
    Dhelp(10,5) = 1.0_DP/8.0_dp*dy
    
    ! WARNING: NOT TESTED!!!

    ! xx-derivatives on current element
    if (Bder(DER_DERIV2D_XX)) then
      dc1 = dxjs * Djac(4)**2
      dc2 = dxjs * ( -2.0_DP * Djac(4) * Djac(2) )
      dc3 = dxjs * Djac(2)**2
      do idof = 1,10
        Dbas(idof,DER_DERIV2D_XX) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif

    ! xy-derivatives on current element
    if (Bder(DER_DERIV2D_XY)) then
      dc1 = - dxjs * Djac(4) * Djac(3)
      dc2 = dxjs * ( Djac(4) * Djac(1) + Djac(3) * Djac(2) )
      dc3 = - dxjs * Djac(1) * Djac(4)
      do idof = 1,10
        Dbas(idof,DER_DERIV2D_XY) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif

    ! yy-derivatives on current element
    if (Bder(DER_DERIV2D_YY)) then
      dc1 = dxjs * Djac(3)**2
      dc2 = dxjs * ( -2.0_DP * Djac(3) * Djac(1) )
      dc3 = dxjs * Djac(1)**2
      do idof = 1,10
        Dbas(idof,DER_DERIV2D_YY) = &
            dc1 * Dhelp(idof,3) + dc2 * Dhelp(idof,4) + dc3 * Dhelp(idof,5)
      end do
    endif
  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T3_2D_mult (celement, Dcoords, Djac, Ddetj, &
                                      Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T3_2D.
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
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder

  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
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

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(10,5,npoints) :: Dhelp
  real(DP) :: dx,dy
  real(DP),dimension(npoints) :: Dxj,Dxjs,Dc1,Dc2,Dc3 !auxiliary variable

  integer :: i,idof   ! point counter

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q3 = 1.0_DP/3.0_DP
  real(DP), parameter :: Q4 = 0.25_DP
  real(DP), parameter :: Q6 = 1.0_DP/6.0_DP
  real(DP), parameter :: Q48 = 1.0_DP/48.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T1_2D-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dbas(1,DER_FUNC,i)= 1.0_DP
      Dbas(2,DER_FUNC,i)= dx
      Dbas(3,DER_FUNC,i)= dy
      Dbas(4,DER_FUNC,i)= Q2*dx*dx-Q6
      Dbas(5,DER_FUNC,i)= Q2*dy*dy-Q6
      Dbas(6,DER_FUNC,i)= dx*dy
      Dbas(7,DER_FUNC,i)= Q48*dx*dx*dx
      Dbas(8,DER_FUNC,i)= Q48*dx*dx*dy
      Dbas(9,DER_FUNC,i)= Q48*dx*dy*dy
      Dbas(10,DER_FUNC,i)= Q48*dy*dy*dy
    end do
  !ENDIF

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    Dxj(:) = 1.0_DP / Ddetj(1:npoints)

    !x- and y-derivatives on reference element
    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)

      !x- and y-derivatives on reference element
      Dhelp(1,1,i)= 0.0_DP
      Dhelp(2,1,i)= 1.0_DP
      Dhelp(3,1,i)= 0.0_DP
      Dhelp(4,1,i)= dx
      Dhelp(5,1,i)= 0.0_DP
      Dhelp(6,1,i)= dy
      Dhelp(7,1,i)= 1.0_dp/16.0_dp*dx*dx
      Dhelp(8,1,i)= 1.0_dp/24.0_dp*dx*dy
      Dhelp(9,1,i)= 1.0_dp/48.0_dp*dy*dy
      Dhelp(10,1,i)= 0.0_dp

      Dhelp(1,2,i)= 0.0_DP
      Dhelp(2,2,i)= 0.0_DP
      Dhelp(3,2,i)= 1.0_DP
      Dhelp(4,2,i)= 0.0_DP
      Dhelp(5,2,i)= dy
      Dhelp(6,2,i)= dx
      Dhelp(7,2,i)= 0.0_dp
      Dhelp(8,2,i)= 1.0_dp/48.0_dp*dx*dx
      Dhelp(9,2,i)= 1.0_dp/24.0_dp*dx*dy
      Dhelp(10,2,i)= 1.0_dp/16.0_dp*dy*dy
      
    end do

    ! x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        do idof = 1,10
          Dbas(idof,DER_DERIV_X,i) = &
              Dxj(i) * (Djac(4,i) * Dhelp(idof,1,i) - Djac(2,i) * Dhelp(idof,2,i))
        end do
!      END DO
!    ENDIF

    ! y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        do idof = 1,10
          Dbas(idof,DER_DERIV_Y,i) = &
              -Dxj(i) * (Djac(3,i) * Dhelp(idof,1,i) - Djac(1,i) * Dhelp(idof,2,i))
        end do
      end do
!    ENDIF
  endif

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]

    !x- and y-derivatives on reference element
    Dxj(:) = 1.0E0_DP / Ddetj(1:npoints)
    Dxjs = Dxj*Dxj

    do i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)

      !xx-, xy and yy-derivatives on reference element
          Dhelp(1,3,i) = 0.0_DP
          Dhelp(2,3,i) = 0.0_DP
          Dhelp(3,3,i) = 0.0_DP
          Dhelp(4,3,i) = 1.0_DP
          Dhelp(5,3,i) = 0.0_DP
          Dhelp(6,3,i) = 0.0_DP
          Dhelp(7,3,i) = 1.0_DP/8.0_dp*dx
          Dhelp(8,3,i) = 1.0_DP/24.0_dp*dy
          Dhelp(9,3,i) = 0.0_DP
          Dhelp(10,3,i) = 0.0_DP

          Dhelp(1,4,i) = 0.0_DP
          Dhelp(2,4,i) = 0.0_DP
          Dhelp(3,4,i) = 0.0_DP
          Dhelp(4,4,i) = 0.0_DP
          Dhelp(5,4,i) = 0.0_DP
          Dhelp(6,4,i) = 1.0_DP
          Dhelp(7,4,i) = 0.0_dp
          Dhelp(8,4,i) = 1.0_DP/24.0_dp*dx
          Dhelp(9,4,i) = 1.0_DP/24.0_dp*dy
          Dhelp(10,4,i) = 0.0_DP

          Dhelp(1,5,i) = 0.0_DP
          Dhelp(2,5,i) = 0.0_DP
          Dhelp(3,5,i) = 0.0_DP
          Dhelp(4,5,i) = 0.0_DP
          Dhelp(5,5,i) = 1.0_DP
          Dhelp(6,5,i) = 0.0_DP
          Dhelp(7,5,i) = 0.0_dp
          Dhelp(8,5,i) = 0.0_dp
          Dhelp(9,5,i) = 1.0_DP/24.0_dp*dx
          Dhelp(10,5,i) = 1.0_DP/8.0_dp*dy

    end do

    ! WARNING: NOT TESTED!!!

    ! xx-derivatives on current element
    !if (Bder(DER_DERIV_XX)) then
      do i=1,npoints
        Dc1(i) = Dxjs(i) * Djac(4,i)**2
        Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(4,i) * Djac(2,i) )
        Dc3(i) = Dxjs(i) * Djac(2,i)**2
      end do

      do i=1,npoints
        do idof = 1,10
          Dbas(idof,DER_DERIV_XX,i) = &
              Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
        end do
      end do
    !endif

    ! xy-derivatives on current element
    !if (Bder(DER_DERIV_XY)) then
      do i=1,npoints
        Dc1(i) = - Dxjs(i) * Djac(4,i) * Djac(3,i)
        Dc2(i) = Dxjs(i) * ( Djac(4,i) * Djac(1,i) + Djac(3,i) * Djac(2,i) )
        Dc3(i) = - Dxjs(i) * Djac(1,i) * Djac(4,i)
      end do

      do i=1,npoints
        do idof = 1,10
          Dbas(idof,DER_DERIV_XY,i) = &
              Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
        end do
      end do
    !endif

    ! yy-derivatives on current element
    !if (Bder(DER_DERIV_YY)) then
      do i=1,npoints
        Dc1(i) = Dxjs(i) * Djac(3,i)**2
        Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(3,i) * Djac(1,i) )
        Dc3(i) = Dxjs(i) * Djac(1,i)**2
      end do

      do idof = 1,10
        Dbas(idof,DER_DERIV_YY,i) = &
            Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
      end do
    !endif
  endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_DG_T3_2D_sim (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, nelements, &
                                Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T3_2D.
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
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Local performance configuration.
  type(t_perfconfig), intent(in) :: rperfconfig
!</input>

!<output>
  ! Value/derivatives of basis functions.
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i-th
  !   basis function of the finite element k in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(10,5,npoints) :: Dhelp
  real(DP) :: dx,dy
  real(DP),dimension(npoints) :: Dxj,Dxjs,Dc1,Dc2,Dc3 !auxiliary variables

  integer :: i   ! point counter
  integer :: j   ! element counter
  integer :: idof

  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q3 = 1.0_DP/3.0_DP
  real(DP), parameter :: Q4 = 0.25_DP
  real(DP), parameter :: Q6 = 1.0_DP/6.0_DP
  real(DP), parameter :: Q48 = 1.0_DP/48.0_DP

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T3_2D-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i,dx,dy) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dbas(1,DER_FUNC,i,j)= 1.0_DP
        Dbas(2,DER_FUNC,i,j)= dx
        Dbas(3,DER_FUNC,i,j)= dy
        Dbas(4,DER_FUNC,i,j)= Q2*dx*dx-Q6
        Dbas(5,DER_FUNC,i,j)= Q2*dy*dy-Q6
        Dbas(6,DER_FUNC,i,j)= dx*dy
        Dbas(7,DER_FUNC,i,j)= Q48*dx*dx*dx
        Dbas(8,DER_FUNC,i,j)= Q48*dx*dx*dy
        Dbas(9,DER_FUNC,i,j)= Q48*dx*dy*dy
        Dbas(10,DER_FUNC,i,j)= Q48*dy*dy*dy
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,Dxj,dx,dy,Dhelp,idof) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1.0E0_DP / Ddetj(1:npoints,j)

      !x- and y-derivatives on reference element
      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)

        !x- and y-derivatives on reference element
        Dhelp(1,1,i)= 0.0_DP
        Dhelp(2,1,i)= 1.0_DP
        Dhelp(3,1,i)= 0.0_DP
        Dhelp(4,1,i)= dx
        Dhelp(5,1,i)= 0.0_DP
        Dhelp(6,1,i)= dy
        Dhelp(7,1,i)= 1.0_dp/16.0_dp*dx*dx
        Dhelp(8,1,i)= 1.0_dp/24.0_dp*dx*dy
        Dhelp(9,1,i)= 1.0_dp/48.0_dp*dy*dy
        Dhelp(10,1,i)= 0.0_dp

        Dhelp(1,2,i)= 0.0_DP
        Dhelp(2,2,i)= 0.0_DP
        Dhelp(3,2,i)= 1.0_DP
        Dhelp(4,2,i)= 0.0_DP
        Dhelp(5,2,i)= dy
        Dhelp(6,2,i)= dx
        Dhelp(7,2,i)= 0.0_dp
        Dhelp(8,2,i)= 1.0_dp/48.0_dp*dx*dx
        Dhelp(9,2,i)= 1.0_dp/24.0_dp*dx*dy
        Dhelp(10,2,i)= 1.0_dp/16.0_dp*dy*dy

      end do

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          do idof = 1,10
            Dbas(idof,DER_DERIV_X,i,j) = &
                Dxj(i) * (Djac(4,i,j) * Dhelp(idof,1,i) &
                          - Djac(2,i,j) * Dhelp(idof,2,i))
!          end do
!        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
!        do i=1,npoints
!          do idof = 1,6
            Dbas(idof,DER_DERIV_Y,i,j) = &
                -Dxj(i) * (Djac(3,i,j) * Dhelp(idof,1,i) &
                - Djac(1,i,j) * Dhelp(idof,2,i))
          end do
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  ! if xx-, xy or yy-derivatives are desired
  if (Bder(DER_DERIV_XX) .or. Bder(DER_DERIV_XY) .or. Bder(DER_DERIV_YY)) then

    !x- and y-derivatives on reference element

    ! Now we have to compute the derivatives on the real element by those on
    ! the reference element. This is a little bit tricky!
    !
    ! Let us assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It is inverse s^{-1} maps to the reference element.
    ! We have:
    !
    !   f(x) = g(y) = g(s^{-1}(x))
    !
    ! We want to compute the hessian
    !
    !   Hf(x) = [ f11 f12 ]
    !           [ f21 f22 ]
    !
    ! with fij = di dj f. By the chain rule, we have:
    !
    !  Hf(x)  =  H [ g(s^{-1}(x)) ]
    !         =  D [ D ( g(s^{-1}(x)) ) ]
    !         =  D [ Dg (s^{-1}(x)) Ds^{-1}(x) ]
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)  +  Dg (s^{-1}(x)) Hs^{-1}(x)
    !
    ! Now the big assumption: We assume that out mapping s is bilinear!
    ! From this assumption we derive that Hs^{-1}(x), so the 2nd term
    ! vanishes. We clearly point out that this is of course not true for
    ! isoparametric mappings!
    !
    ! With this assumption, we continue:
    !
    !         =  D [ Dg (s^{-1}(x)) ] Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(s^{-1}(x))  Ds^{-1}(x)
    !         =  Ds^{-1}(x)^T  Hg(y)  Ds^{-1}(x)
    !
    ! This is computable now. Let us write:
    !
    !  Ds = [ s11 s12 ] , Ds^{-1} = 1/det [  s22 -s12 ] , Hg(y) = [ g11 g12 ]
    !       [ s21 s22 ]                   [ -s21  s11 ]           [ g21 g22 ]
    !
    ! then we get:
    !
    !  Hf(x) = 1/det^2
    !  [(s22*g11-s21*g21)*s22-(s22*g12-s21*g22)*s21, -(s22*g11-s21*g21)*s12+(s22*g12-s21*g22)*s11],
    !  [(-s12*g11+s11*g21)*s22-(-s12*g12+s11*g22)*s21, -(-s12*g11+s11*g21)*s12+(-s12*g12+s11*g22)*s11]])
    !
    ! This can be simplified by the fact that g12=g21 (as (gij) is a Hessian):
    !
    !  = 1/det^2
    !    [s22^2*g11-2*s22*s21*g21+s21^2*g22, -s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22]
    !    [-s22*s12*g11+s22*s11*g21+s12*s21*g21-s21*s11*g22, s12^2*g11-2*s12*s11*g21+s11^2*g22]])
    !
    ! so we have
    !
    !  f11       = 1/det^2 * [ s22^2*g11    - 2*s22*s21*g21               + s21^2*g22 ]
    !  f12 = f21 = 1/det^2 * [ -s22*s12*g11 + ( s22*s11 + s12*s21 ) * g21 - s21*s11*g22 ]
    !  f22       = 1/det^2 * [ s12^2*g11    - 2*s12*s11*g21               + s11^2*g22 ]

    !x- and y-derivatives on reference element
    !$omp parallel do default(shared)&
    !$omp private(Dc1,Dc2,Dc3,Dhelp,Dxj,Dxjs,dx,dy,i,idof)&
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1.0E0_DP / Ddetj(1:npoints,j)
      Dxjs = Dxj*Dxj

      do i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)

        !xx-, xy and yy-derivatives on reference element
        Dhelp(1,3,i) = 0.0_DP
        Dhelp(2,3,i) = 0.0_DP
        Dhelp(3,3,i) = 0.0_DP
        Dhelp(4,3,i) = 1.0_DP
        Dhelp(5,3,i) = 0.0_DP
        Dhelp(6,3,i) = 0.0_DP
        Dhelp(7,3,i) = 1.0_DP/8.0_dp*dx
        Dhelp(8,3,i) = 1.0_DP/24.0_dp*dy
        Dhelp(9,3,i) = 0.0_DP
        Dhelp(10,3,i) = 0.0_DP

        Dhelp(1,4,i) = 0.0_DP
        Dhelp(2,4,i) = 0.0_DP
        Dhelp(3,4,i) = 0.0_DP
        Dhelp(4,4,i) = 0.0_DP
        Dhelp(5,4,i) = 0.0_DP
        Dhelp(6,4,i) = 1.0_DP
        Dhelp(7,4,i) = 0.0_dp
        Dhelp(8,4,i) = 1.0_DP/24.0_dp*dx
        Dhelp(9,4,i) = 1.0_DP/24.0_dp*dy
        Dhelp(10,4,i) = 0.0_DP

        Dhelp(1,5,i) = 0.0_DP
        Dhelp(2,5,i) = 0.0_DP
        Dhelp(3,5,i) = 0.0_DP
        Dhelp(4,5,i) = 0.0_DP
        Dhelp(5,5,i) = 1.0_DP
        Dhelp(6,5,i) = 0.0_DP
        Dhelp(7,5,i) = 0.0_dp
        Dhelp(8,5,i) = 0.0_dp
        Dhelp(9,5,i) = 1.0_DP/24.0_dp*dx
        Dhelp(10,5,i) = 1.0_DP/8.0_dp*dy

      end do

      ! WARNING: NOT TESTED!!!

      ! xx-derivatives on current element
      !if (Bder(DER_DERIV_XX)) then
        do i=1,npoints
          Dc1(i) = Dxjs(i) * Djac(4,i,j)**2
          Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(4,i,j) * Djac(2,i,j) )
          Dc3(i) = Dxjs(i) * Djac(2,i,j)**2
        end do

        do i=1,npoints
          do idof = 1,10
            Dbas(idof,DER_DERIV_XX,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

      ! xy-derivatives on current element
      !if (Bder(DER_DERIV_XY)) then
        do i=1,npoints
          Dc1(i) = - Dxjs(i) * Djac(4,i,j) * Djac(3,i,j)
          Dc2(i) = Dxjs(i) * ( Djac(4,i,j) * Djac(1,i,j) + Djac(3,i,j) * Djac(2,i,j) )
          Dc3(i) = - Dxjs(i) * Djac(1,i,j) * Djac(4,i,j)
        end do

        do i=1,npoints
          do idof = 1,10
            Dbas(idof,DER_DERIV_XY,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

      ! yy-derivatives on current element
      !if (Bder(DER_DERIV_YY)) then
        do i=1,npoints
          Dc1(i) = Dxjs(i) * Djac(3,i,j)**2
          Dc2(i) = Dxjs(i) * ( -2.0_DP * Djac(3,i,j) * Djac(1,i,j) )
          Dc3(i) = Dxjs(i) * Djac(1,i,j)**2
        end do

        do i=1,npoints
          do idof = 1,10
            Dbas(idof,DER_DERIV_YY,i,j) = &
                Dc1(i) * Dhelp(idof,3,i) + Dc2(i) * Dhelp(idof,4,i) + Dc3(i) * Dhelp(idof,5,i)
          end do
        end do
      !endif

    end do
    !$omp end parallel do

  end if

  end subroutine

  !****************************************************************************
  !****************************************************************************

  ! -------------- NEW ELEMENT INTERFACE IMPLEMENTATIONS FOLLOW --------------

  !****************************************************************************
  !****************************************************************************



  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_Q1_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The Q1_2D element is specified by four polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,4:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Pi(vj) = kronecker(i,j)
  !   }
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the quadrilateral
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  ! P1(x,y) = 1/4 * (1 - x) * (1 - y)
  ! P2(x,y) = 1/4 * (1 + x) * (1 - y)
  ! P3(x,y) = 1/4 * (1 + x) * (1 + y)
  ! P4(x,y) = 1/4 * (1 - x) * (1 + y)

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 4

  ! Local variables
  real(DP) :: ddet,dx,dy
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! If function values are desired, calculate them.
      !
      ! We have to compute the basis functions in the points.
      ! I.e., we have to evaluate
      !
      !    phi_k(x) = Pk(sigma^-1(x)) = Pk^(x^)
      !
      ! with x being the real world coordinates, x^ the coordinates
      ! in the reference element and sigma: x^ -> x the mapping
      ! between the reference and the real element.
      ! So just evaluate the basis functions in the points x^ 
      ! on the reference element

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP-dx)*(1.0_DP-dy)
          Dbas(2,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP+dx)*(1.0_DP-dy)
          Dbas(3,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP+dx)*(1.0_DP+dy)
          Dbas(4,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP-dx)*(1.0_DP+dy)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! If x-or y-derivatives are desired, calculate them.
      ! The values of the derivatives are calculated by taking the
      ! derivative of the polynomials and multiplying them with the
      ! inverse of the transformation matrix (in each point) as
      ! stated above.
      !
      ! We have to evaluate "grad(phi(x))". This is done by
      ! using the chain rule as follows:
      !
      !    grad(phi_k(x))^T = D phi_k(x)
      !                     = D Pk(sigma^-1(x))
      !                     = D Pk(x^) * D sigma^-1(x)
      !
      ! Now note that the Jacobian "D sigma(x^)" of the mapping
      ! between the reference and the real element is given in Djac:
      !
      !    D sigma(x^) = ( Djac(1) Djac(3) )
      !                  ( Djac(2) Djac(4) )
      !
      ! Its inverse then reads
      !
      !    D sigma^-1(x) = 1/det ( -Djac(4)  Djac(3) )
      !                          (  Djac(2) -Djac(1) )
      !
      ! with det = Djac(1)*Djac(4) - Djac(2)*Djac(3).
      ! So all in all, the derivative can be computed as follows:
      !
      !    grad(phi_k(x))^T = D Pk(x^) * D sigma^-1(x)
      !                     = ( dx(Pk(x^)) dy(Pk(x^)) ) * 1/det ( -Djac(4)  Djac(3) )
      !                                                         (  Djac(2) -Djac(1) )
      ! i.e., we have
      !
      !   dx(phi_k(x)) = 1/det ( -dx(Pk(x^))*Djac(4) + dy(Pk(x^))*Djac(2) )
      !   dy(phi_k(x)) = 1/det (  dx(Pk(x^))*Djac(3) - dy(Pk(x^))*Djac(1) )

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) = -0.25_DP*(1.0_DP-dy)
          DrefDer(2,1) =  0.25_DP*(1.0_DP-dy)
          DrefDer(3,1) =  0.25_DP*(1.0_DP+dy)
          DrefDer(4,1) = -0.25_DP*(1.0_DP+dy)
          ! Y-derivatives
          DrefDer(1,2) = -0.25_DP*(1.0_DP-dx)
          DrefDer(2,2) = -0.25_DP*(1.0_DP+dx)
          DrefDer(3,2) =  0.25_DP*(1.0_DP+dx)
          DrefDer(4,2) =  0.25_DP*(1.0_DP-dx)

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_Q1B_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The Q1B_2D element is specified by five polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y , (1-x^2)*(1-y^2) }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,4:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Pi(vj) = kronecker(i,j)
  !   }
  !   Pi(x) = kronecker(i,5)
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the quadrilateral
  ! x beind the midpoint of the quadrilateral
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  ! P1(x,y) = 1/4 * ((1 - x) * (1 - y) - (1 - x^2) * (1 - y^2))
  ! P2(x,y) = 1/4 * ((1 + x) * (1 - y) - (1 - x^2) * (1 - y^2))
  ! P3(x,y) = 1/4 * ((1 + x) * (1 + y) - (1 - x^2) * (1 - y^2))
  ! P4(x,y) = 1/4 * ((1 - x) * (1 + y) - (1 - x^2) * (1 - y^2))
  ! P5(x,y) = (1 - x^2) * (1 - y^2)

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 5

  ! Local variables
  real(DP) :: ddet,dx,dy,db
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! If function values are desired, calculate them.
      !
      ! We have to compute the basis functions in the points.
      ! I.e., we have to evaluate
      !
      !    phi_k(x) = Pk(sigma^-1(x)) = Pk^(x^)
      !
      ! with x being the real world coordinates, x^ the coordinates
      ! in the reference element and sigma: x^ -> x the mapping
      ! between the reference and the real element.
      ! So just evaluate the basis functions in the points x^ 
      ! on the reference element

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)
          
          ! compute bubble
          db = (1.0_DP - dx**2) * (1.0_DP - dy**2)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = 0.25_DP*((1.0_DP-dx)*(1.0_DP-dy)-db)
          Dbas(2,DER_FUNC2D,i,j) = 0.25_DP*((1.0_DP+dx)*(1.0_DP-dy)-db)
          Dbas(3,DER_FUNC2D,i,j) = 0.25_DP*((1.0_DP+dx)*(1.0_DP+dy)-db)
          Dbas(4,DER_FUNC2D,i,j) = 0.25_DP*((1.0_DP-dx)*(1.0_DP+dy)-db)
          Dbas(5,DER_FUNC2D,i,j) = db

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! If x-or y-derivatives are desired, calculate them.
      ! The values of the derivatives are calculated by taking the
      ! derivative of the polynomials and multiplying them with the
      ! inverse of the transformation matrix (in each point) as
      ! stated above.
      !
      ! We have to evaluate "grad(phi(x))". This is done by
      ! using the chain rule as follows:
      !
      !    grad(phi_k(x))^T = D phi_k(x)
      !                     = D Pk(sigma^-1(x))
      !                     = D Pk(x^) * D sigma^-1(x)
      !
      ! Now note that the Jacobian "D sigma(x^)" of the mapping
      ! between the reference and the real element is given in Djac:
      !
      !    D sigma(x^) = ( Djac(1) Djac(3) )
      !                  ( Djac(2) Djac(4) )
      !
      ! Its inverse then reads
      !
      !    D sigma^-1(x) = 1/det ( -Djac(4)  Djac(3) )
      !                          (  Djac(2) -Djac(1) )
      !
      ! with det = Djac(1)*Djac(4) - Djac(2)*Djac(3).
      ! So all in all, the derivative can be computed as follows:
      !
      !    grad(phi_k(x))^T = D Pk(x^) * D sigma^-1(x)
      !                     = ( dx(Pk(x^)) dy(Pk(x^)) ) * 1/det ( -Djac(4)  Djac(3) )
      !                                                         (  Djac(2) -Djac(1) )
      ! i.e., we have
      !
      !   dx(phi_k(x)) = 1/det ( -dx(Pk(x^))*Djac(4) + dy(Pk(x^))*Djac(2) )
      !   dy(phi_k(x)) = 1/det (  dx(Pk(x^))*Djac(3) - dy(Pk(x^))*Djac(1) )

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) = -0.25_DP*(1.0_DP-dy)+0.5_DP*dx*(1.0_DP-dy**2)
          DrefDer(2,1) =  0.25_DP*(1.0_DP-dy)+0.5_DP*dx*(1.0_DP-dy**2)
          DrefDer(3,1) =  0.25_DP*(1.0_DP+dy)+0.5_DP*dx*(1.0_DP-dy**2)
          DrefDer(4,1) = -0.25_DP*(1.0_DP+dy)+0.5_DP*dx*(1.0_DP-dy**2)
          DrefDer(5,1) = -2.0_DP*dx*(1.0_DP-dy**2)
          ! Y-derivatives
          DrefDer(1,2) = -0.25_DP*(1.0_DP-dx)+0.5_DP*dy*(1.0_DP-dx**2)
          DrefDer(2,2) = -0.25_DP*(1.0_DP+dx)+0.5_DP*dy*(1.0_DP-dx**2)
          DrefDer(3,2) =  0.25_DP*(1.0_DP+dx)+0.5_DP*dy*(1.0_DP-dx**2)
          DrefDer(4,2) =  0.25_DP*(1.0_DP-dx)+0.5_DP*dy*(1.0_DP-dx**2)
          DrefDer(5,2) = -2.0_DP*dy*(1.0_DP-dx**2)

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_EM11_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The EM11_2D element is specified by four polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,4:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Pi(vj) = kronecker(i,j)
  !   }
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the quadrilateral


  ! Parameter: Number of local basis functions
  integer, parameter :: NBAS = 4

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr
  real(DP) :: ddets

  ! other local variables
  logical :: bsuccess
  integer :: i,iel,ipt
  real(DP), dimension(NBAS,NBAS) :: Da, Dc
  real(DP) :: dx,dy,derx,dery

    ! Loop over all elements
    !$omp parallel do default(shared) &
    !$omp private(Da,Dc,Dr,Ds,Dvert,bsuccess,ddets,derx,dery,dx,dy,i,ipt) &
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Fetch vertice coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - dr(2)
      ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - dr(1))
      ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - dr(2))
      ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - dr(1)
      ddets = 1.0_DP / (ds(1,1)*ds(2,2) - ds(1,2)*ds(2,1))
      Ds = ddets * Ds

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Build coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Clear coefficient matrix
      Da = 0.0_DP

      ! Loop over all vertices of the quad
      do i = 1, 4

        ! Apply inverse affine trafo on corner vertice coordinates to get
        ! (x,y) in the element`s local coordinate system.
        dx = ds(1,1)*(Dvert(1,i)-dr(1)) + ds(1,2)*(Dvert(2,i)-dr(2))
        dy = ds(2,1)*(Dvert(1,i)-dr(1)) + ds(2,2)*(Dvert(2,i)-dr(2))

        ! Evaluate monomials in current vertice
        Da(1,i) = 1.0_DP
        Da(2,i) = dx
        Da(3,i) = dy
        Da(4,i) = dx*dy

      end do ! i

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Invert coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Call the 'direct' inversion routine for 4x4 systems
      call mprim_invert4x4MatrixDirect(Da, Dc, bsuccess)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 5: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate basis functions
          do i = 1, NBAS

            Dbas(i,DER_FUNC2D,ipt,iel) = Dc(i,1) + dx*Dc(i,2) &
                                       + dy*(Dc(i,3) + dx*Dc(i,4))
          end do ! i

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 6: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate derivatives
          do i = 1, NBAS

            ! Calculate 'reference' derivatives
            derx = Dc(i,2) + dy*Dc(i,4)
            dery = Dc(i,3) + dx*Dc(i,4)

            ! Calculate 'real' derivatives
            Dbas(i,DER_DERIV2D_X,ipt,iel) = ds(1,1)*derx + ds(2,1)*dery
            Dbas(i,DER_DERIV2D_Y,ipt,iel) = ds(1,2)*derx + ds(2,2)*dery

          end do ! i

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_Q2_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The Q2_2D element is specified by nine polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2, y^2, x^2*y, x*y^2, x^2*y^2}
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,9:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Pi(vj) = kronecker(i,j)
  !     Pi(ej) = kronecker(i,j+4)
  !   }
  !   Pi(0,0) = kronecker(i,9)
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the quadrilateral
  ! ej being the midpoint of the j-th local edge of the quadrilateral
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  ! P1(x,y) =  1/4 * (1 - x) * (1 - y) * x * y
  ! P2(x,y) = -1/4 * (1 + x) * (1 - y) * x * y
  ! P3(x,y) =  1/4 * (1 + x) * (1 + y) * x * y
  ! P4(x,y) = -1/4 * (1 - x) * (1 + y) * x * y
  ! P5(x,y) = -1/2 * (1 - x^2) * (1 - y  ) * y
  ! P6(x,y) =  1/2 * (1 + x  ) * (1 - y^2) * x
  ! P7(x,y) =  1/2 * (1 - x^2) * (1 + y  ) * y
  ! P8(x,y) = -1/2 * (1 - x  ) * (1 - y^2) * x
  ! P9(x,y) =  (1 - x^2) * (1 - y^2)

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 9

  ! Local variables
  real(DP) :: ddet,dx,dy
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) =  0.25_DP*(1.0_DP-dx)*(1.0_DP-dy)*dx*dy
          Dbas(2,DER_FUNC2D,i,j) = -0.25_DP*(1.0_DP+dx)*(1.0_DP-dy)*dx*dy
          Dbas(3,DER_FUNC2D,i,j) =  0.25_DP*(1.0_DP+dx)*(1.0_DP+dy)*dx*dy
          Dbas(4,DER_FUNC2D,i,j) = -0.25_DP*(1.0_DP-dx)*(1.0_DP+dy)*dx*dy
          Dbas(5,DER_FUNC2D,i,j) = -0.5_DP*(1.0_DP-dx*dx)*(1.0_DP-dy)*dy
          Dbas(6,DER_FUNC2D,i,j) =  0.5_DP*(1.0_DP+dx)*(1.0_DP-dy*dy)*dx
          Dbas(7,DER_FUNC2D,i,j) =  0.5_DP*(1.0_DP-dx*dx)*(1.0_DP+dy)*dy
          Dbas(8,DER_FUNC2D,i,j) = -0.5_DP*(1.0_DP-dx)*(1.0_DP-dy*dy)*dx
          Dbas(9,DER_FUNC2D,i,j) = (1.0_DP-dx*dx)*(1.0_DP-dy*dy)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) =  0.25_DP*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy)*dy
          DrefDer(2,1) = -0.25_DP*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy)*dy
          DrefDer(3,1) =  0.25_DP*(1.0_DP+2.0_DP*dx)*(1.0_DP+dy)*dy
          DrefDer(4,1) = -0.25_DP*(1.0_DP-2.0_DP*dx)*(1.0_DP+dy)*dy
          DrefDer(5,1) =  (1.0_DP-dy)*dx*dy
          DrefDer(6,1) =  0.5_DP*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy*dy)
          DrefDer(7,1) = -(1.0_DP+dy)*dx*dy
          DrefDer(8,1) = -0.5_DP*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy*dy)
          DrefDer(9,1) = -2.0_DP*(1.0_DP-dy*dy)*dx
          ! Y-derivatives
          DrefDer(1,2) =  0.25_DP*(1.0_DP-dx)*(1.0_DP-2.0_DP*dy)*dx
          DrefDer(2,2) = -0.25_DP*(1.0_DP+dx)*(1.0_DP-2.0_DP*dy)*dx
          DrefDer(3,2) =  0.25_DP*(1.0_DP+dx)*(1.0_DP+2.0_DP*dy)*dx
          DrefDer(4,2) = -0.25_DP*(1.0_DP-dx)*(1.0_DP+2.0_DP*dy)*dx
          DrefDer(5,2) = -0.5_DP*(1.0_DP-dx*dx)*(1.0_DP-2.0_DP*dy)
          DrefDer(6,2) = -(1.0_DP+dx)*dx*dy
          DrefDer(7,2) =  0.5_DP*(1.0_DP-dx*dx)*(1.0_DP+2.0_DP*dy)
          DrefDer(8,2) =  (1.0_DP-dx)*dx*dy
          DrefDer(9,2) = -2.0_DP*(1.0_DP-dx*dx)*dy

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_Q2H_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! THIS ELEMENT IS EXPERIMENTAL !!!

  ! Element Description
  ! -------------------
  ! The Q2H_2D element is specified by nine polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2, y^2, x^2*y, x*y^2, x^2*y^2}
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,9:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Pi(vj) = kronecker(i,j)
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))*L2(t)) d(t) = kronecker(i,j+4) * |ej|
  !   }
  !   Int_T (Pi(x,y)*L2(x)*L2(y)) d(x,y) = kronecker(i,9) * |T|
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the quadrilateral
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  ! T being the quadrilateral
  ! |T| being the area of the quadrilateral
  ! L2 being the second Legendre-Polynomial:
  ! L2(x) := 1/2*(3*x^2 - 1)
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  ! P1(x,y) = 1/4 * (1 - x) * (1 - y)
  ! P2(x,y) = 1/4 * (1 + x) * (1 - y)
  ! P3(x,y) = 1/4 * (1 + x) * (1 + y)
  ! P4(x,y) = 1/4 * (1 - x) * (1 + y)
  ! P5(x,y) = 15/8 * (x^2 * (1 - y) + y - 1)
  ! P6(x,y) = 15/8 * (y^2 * (1 + x) - x - 1)
  ! P7(x,y) = 15/8 * (x^2 * (1 + y) - y - 1)
  ! P8(x,y) = 15/8 * (y^2 * (1 - x) + x - 1)
  ! P9(x,y) = -225/16 * (x^2 + y^2 - x^2*y^2 - 1)

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 9

  ! Local variables
  real(DP) :: ddet,dx,dy,dx2,dy2
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy,dx2,dy2)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Pre-compute squares
          dx2 = dx*dx
          dy2 = dy*dy

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP-dx)*(1.0_DP-dy)
          Dbas(2,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP+dx)*(1.0_DP-dy)
          Dbas(3,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP+dx)*(1.0_DP+dy)
          Dbas(4,DER_FUNC2D,i,j) = 0.25_DP*(1.0_DP-dx)*(1.0_DP+dy)
          Dbas(5,DER_FUNC2D,i,j) = 1.875_DP*(dx2*(1.0_DP - dy) + dy - 1.0_DP)
          Dbas(6,DER_FUNC2D,i,j) = 1.875_DP*(dy2*(1.0_DP + dx) - dx - 1.0_DP)
          Dbas(7,DER_FUNC2D,i,j) = 1.875_DP*(dx2*(1.0_DP + dy) - dy - 1.0_DP)
          Dbas(8,DER_FUNC2D,i,j) = 1.875_DP*(dy2*(1.0_DP - dx) + dx - 1.0_DP)
          Dbas(9,DER_FUNC2D,i,j) = -14.0625 * (dx2 + dy2 - dx2*dy2 - 1.0_DP)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy,dx2,dy2)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Pre-compute squares
          dx2 = dx*dx
          dy2 = dy*dy

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) = -0.25_DP*(1.0_DP-dy)
          DrefDer(2,1) =  0.25_DP*(1.0_DP-dy)
          DrefDer(3,1) =  0.25_DP*(1.0_DP+dy)
          DrefDer(4,1) = -0.25_DP*(1.0_DP+dy)
          DrefDer(5,1) =  3.75_DP*dx*(1.0_DP - dy)
          DrefDer(6,1) = -1.875_DP*(1.0_DP - dy2)
          DrefDer(7,1) =  3.75_DP*dx*(1.0_DP + dy)
          DrefDer(8,1) =  1.875_DP*(1.0_DP - dy2)
          DrefDer(9,1) = -28.125_DP*dx*(1 - dy2)
          ! Y-derivatives
          DrefDer(1,2) = -0.25_DP*(1.0_DP-dx)
          DrefDer(2,2) = -0.25_DP*(1.0_DP+dx)
          DrefDer(3,2) =  0.25_DP*(1.0_DP+dx)
          DrefDer(4,2) =  0.25_DP*(1.0_DP-dx)
          DrefDer(5,2) =  1.875_DP*(1.0_DP - dx2)
          DrefDer(6,2) =  3.75_DP*dy*(1.0_DP + dx)
          DrefDer(7,2) = -1.875_DP*(1.0_DP - dx2)
          DrefDer(8,2) =  3.75_DP*dy*(1.0_DP - dx)
          DrefDer(9,2) = -28.125_DP*dy*(1 - dx2)

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_Q3_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The Q3_2D element is specified by sixteen polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, x^2, x^3 } * { 1, y, y^2, y^3 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,9:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Pi(vj)   = kronecker(i,j)
  !     Pi(ej_1) = kronecker(4+(2*i)-1,j)
  !     Pi(ej_2) = kronecker(4+(2*i)  ,j)
  !     Pi(xj)   = kronecker(12+i, j)
  !   }
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the quadrilateral
  ! ej_1/2 being the first/second Lobatto-Point the j-th local edge of the quadrilateral
  ! x_j being the j-th Lobatto-Point on the quadrilateral
  !
  ! The basis polynomials Pi are constructed as the product
  !
  !    P_{4*i+j} (x,y) := fi(x) * fj(y)
  !
  ! of the 1D Lagrange-Polynomials to the 4 Gauss-Lobatto points on the interval [-1,1],
  ! namely:
  !
  !   f1(x) := 1/8 * (-1 + x*( 1 + x*5*(1 - x)))
  !   f2(x) := 1/8 * (-1 + x*(-1 + x*5*(1 + x)))
  !   f3(x) := 5/8 * (1 + x*(-sqrt(5) + x*(-1 + x*sqrt(5))))
  !   f4(x) := 5/8 * (1 + x*( sqrt(5) + x*(-1 - x*sqrt(5))))

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 16

  ! Local variables
  real(DP) :: ddet,dx,dy
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas( 1,DER_FUNC2D,i,j) = f1(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          Dbas( 2,DER_FUNC2D,i,j) = f2(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          Dbas( 3,DER_FUNC2D,i,j) = f2(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          Dbas( 4,DER_FUNC2D,i,j) = f1(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          Dbas( 5,DER_FUNC2D,i,j) = f3(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          Dbas( 6,DER_FUNC2D,i,j) = f4(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          Dbas( 7,DER_FUNC2D,i,j) = f2(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          Dbas( 8,DER_FUNC2D,i,j) = f2(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))
          Dbas( 9,DER_FUNC2D,i,j) = f4(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          Dbas(10,DER_FUNC2D,i,j) = f3(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          Dbas(11,DER_FUNC2D,i,j) = f1(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))
          Dbas(12,DER_FUNC2D,i,j) = f1(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          Dbas(13,DER_FUNC2D,i,j) = f3(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          Dbas(14,DER_FUNC2D,i,j) = f4(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          Dbas(15,DER_FUNC2D,i,j) = f4(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))
          Dbas(16,DER_FUNC2D,i,j) = f3(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer( 1,1) = df1(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          DrefDer( 2,1) = df2(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          DrefDer( 3,1) = df2(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          DrefDer( 4,1) = df1(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          DrefDer( 5,1) = df3(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          DrefDer( 6,1) = df4(reval%p_DpointsRef(1,i,j)) * f1(reval%p_DpointsRef(2,i,j))
          DrefDer( 7,1) = df2(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          DrefDer( 8,1) = df2(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))
          DrefDer( 9,1) = df4(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          DrefDer(10,1) = df3(reval%p_DpointsRef(1,i,j)) * f2(reval%p_DpointsRef(2,i,j))
          DrefDer(11,1) = df1(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))
          DrefDer(12,1) = df1(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          DrefDer(13,1) = df3(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          DrefDer(14,1) = df4(reval%p_DpointsRef(1,i,j)) * f3(reval%p_DpointsRef(2,i,j))
          DrefDer(15,1) = df4(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))
          DrefDer(16,1) = df3(reval%p_DpointsRef(1,i,j)) * f4(reval%p_DpointsRef(2,i,j))

          ! Y-derivatives
          DrefDer( 1,2) = f1(reval%p_DpointsRef(1,i,j)) * df1(reval%p_DpointsRef(2,i,j))
          DrefDer( 2,2) = f2(reval%p_DpointsRef(1,i,j)) * df1(reval%p_DpointsRef(2,i,j))
          DrefDer( 3,2) = f2(reval%p_DpointsRef(1,i,j)) * df2(reval%p_DpointsRef(2,i,j))
          DrefDer( 4,2) = f1(reval%p_DpointsRef(1,i,j)) * df2(reval%p_DpointsRef(2,i,j))
          DrefDer( 5,2) = f3(reval%p_DpointsRef(1,i,j)) * df1(reval%p_DpointsRef(2,i,j))
          DrefDer( 6,2) = f4(reval%p_DpointsRef(1,i,j)) * df1(reval%p_DpointsRef(2,i,j))
          DrefDer( 7,2) = f2(reval%p_DpointsRef(1,i,j)) * df3(reval%p_DpointsRef(2,i,j))
          DrefDer( 8,2) = f2(reval%p_DpointsRef(1,i,j)) * df4(reval%p_DpointsRef(2,i,j))
          DrefDer( 9,2) = f4(reval%p_DpointsRef(1,i,j)) * df2(reval%p_DpointsRef(2,i,j))
          DrefDer(10,2) = f3(reval%p_DpointsRef(1,i,j)) * df2(reval%p_DpointsRef(2,i,j))
          DrefDer(11,2) = f1(reval%p_DpointsRef(1,i,j)) * df4(reval%p_DpointsRef(2,i,j))
          DrefDer(12,2) = f1(reval%p_DpointsRef(1,i,j)) * df3(reval%p_DpointsRef(2,i,j))
          DrefDer(13,2) = f3(reval%p_DpointsRef(1,i,j)) * df3(reval%p_DpointsRef(2,i,j))
          DrefDer(14,2) = f4(reval%p_DpointsRef(1,i,j)) * df3(reval%p_DpointsRef(2,i,j))
          DrefDer(15,2) = f4(reval%p_DpointsRef(1,i,j)) * df4(reval%p_DpointsRef(2,i,j))
          DrefDer(16,2) = f3(reval%p_DpointsRef(1,i,j)) * df4(reval%p_DpointsRef(2,i,j))

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if
    
    contains
    
      elemental real(DP) function f1(dx)
      real(DP), intent(in) :: dx
        f1 = 0.125_DP * (-1.0_DP + dx*( 1.0_DP + dx*5.0_DP*( 1.0_DP - dx)))
      end function

      elemental real(DP) function f2(dx)
      real(DP), intent(in) :: dx
        f2 = 0.125_DP * (-1.0_DP + dx*(-1.0_DP + dx*5.0_DP*( 1.0_DP + dx)))
      end function

      elemental real(DP) function f3(dx)
      real(DP), intent(in) :: dx
        f3 = 0.625_DP * ( 1.0_DP + dx*(-sqrt(5.0_DP) + dx*(-1.0_DP + dx*sqrt(5.0_DP))))
      end function

      elemental real(DP) function f4(dx)
      real(DP), intent(in) :: dx
        f4 = 0.625_DP * ( 1.0_DP + dx*( sqrt(5.0_DP) + dx*(-1.0_DP - dx*sqrt(5.0_DP))))
      end function
      
      elemental real(DP) function df1(dx)
      real(DP), intent(in) :: dx
        df1 = 0.125_DP * ( 1.0_DP + dx*(10.0_DP - dx*15.0_DP))
      end function

      elemental real(DP) function df2(dx)
      real(DP), intent(in) :: dx
        df2 = 0.125_DP * (-1.0_DP + dx*(10.0_DP + dx*15.0_DP))
      end function

      elemental real(DP) function df3(dx)
      real(DP), intent(in) :: dx
        df3 = -0.625_DP * (sqrt(5.0_DP) + dx*( 2.0_DP - dx*sqrt(45.0_DP)))
      end function

      elemental real(DP) function df4(dx)
      real(DP), intent(in) :: dx
        df4 =  0.625_DP * (sqrt(5.0_DP) + dx*(-2.0_DP - dx*sqrt(45.0_DP)))
      end function
      
    end subroutine
  
  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_QP1_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The QP1_2D element is specified by three polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  ! { 1, x, y }
  !
  ! As the QP1 element is discontinous, the basis polynomials do not have to
  ! fulfill any special conditions - they are simply defined as:
  !
  !  P1 (x,y,z) = 1
  !  P2 (x,y,z) = x
  !  P3 (x,y,z) = y

  ! Local variables
  real(DP) :: ddet
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = 1.0_DP
          Dbas(2,DER_FUNC2D,i,j) = reval%p_DpointsRef(1,i,j)
          Dbas(3,DER_FUNC2D,i,j) = reval%p_DpointsRef(2,i,j)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,ddet)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          Dbas(1,DER_DERIV2D_X,i,j) =  0.0_DP
          Dbas(2,DER_DERIV2D_X,i,j) =  reval%p_Djac(4,i,j)*ddet
          Dbas(3,DER_DERIV2D_X,i,j) = -reval%p_Djac(2,i,j)*ddet

          ! Y-derivatives on real element
          Dbas(1,DER_DERIV2D_Y,i,j) =  0.0_DP
          Dbas(2,DER_DERIV2D_Y,i,j) = -reval%p_Djac(3,i,j)*ddet
          Dbas(3,DER_DERIV2D_Y,i,j) =  reval%p_Djac(1,i,j)*ddet

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_QPW4P0_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! This element is piecewisely defined by 4 triangles in one
  ! quad. There are 5 DOFs: the four corners plus the element
  ! midpoint:
  !
  ! -1,1              1,1       a3                a2
  !   +-------------+             +-------------+
  !   | \    T3   / |             | \         / |
  !   |   \     /   |             |   \     /   |
  !   |     \ /     |             |     \ /     |
  !   | T4   X   T2 |     ->      |      X MP   |
  !   |     / \     |             |     / \     |
  !   |   /     \   |             |   /     \   |
  !   | /    T1   \ |             | /         \ |
  !   +-------------+             +-------------+
  ! -1,-1             1,-1      a0                a1
  !
  ! Therefore, for every point we have to decide in which
  ! sub-triangle we are.
  !
  ! The point ia in one of the subtriangles T1, T2, T3 or T4.
  ! Calculate the line a0->a2 and chech whether it is left
  ! or right. Calculate the line a1->a3 and check whether it is
  ! left or right. That way, we know in which subtriangle the
  ! point is.

  ! Local variables
  real(DP) :: dx,dy
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas(1,DER_FUNC2D,i,j) = 1.0_DP
              Dbas(2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(4,DER_FUNC2D,i,j) = 0.0_DP

            else
              ! We are in T2

              Dbas(1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(2,DER_FUNC2D,i,j) = 1.0_DP
              Dbas(3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(4,DER_FUNC2D,i,j) = 0.0_DP

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas(1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(3,DER_FUNC2D,i,j) = 1.0_DP
              Dbas(4,DER_FUNC2D,i,j) = 0.0_DP

            else
              ! We are in T4

              Dbas(1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(4,DER_FUNC2D,i,j) = 1.0_DP

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_QPW4P1_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! This element is piecewisely defined by 4 triangles in one
  ! quad. There are 5 DOFs: the four corners plus the element
  ! midpoint:
  !
  ! -1,1              1,1       a3                a2
  !   +-------------+             +-------------+
  !   | \    T3   / |             | \         / |
  !   |   \     /   |             |   \     /   |
  !   |     \ /     |             |     \ /     |
  !   | T4   X   T2 |     ->      |      X MP   |
  !   |     / \     |             |     / \     |
  !   |   /     \   |             |   /     \   |
  !   | /    T1   \ |             | /         \ |
  !   +-------------+             +-------------+
  ! -1,-1             1,-1      a0                a1
  !
  ! Therefore, for every point we have to decide in which
  ! sub-triangle we are.
  !
  ! The point ia in one of the subtriangles T1, T2, T3 or T4.
  ! Calculate the line a0->a2 and chech whether it is left
  ! or right. Calculate the line a1->a3 and check whether it is
  ! left or right. That way, we know in which subtriangle the
  ! point is.

  ! Local variables
  real(DP) :: ddet,dx,dy,a11,a12,a21,a22
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas(1,DER_FUNC2D,i,j) = -0.5_DP*dx-0.5_DP*dy
              Dbas(2,DER_FUNC2D,i,j) =  0.5_DP*dx-0.5_DP*dy
              Dbas(3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(4,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(5,DER_FUNC2D,i,j) = 1.0_DP + dy

            else
              ! We are in T2

              Dbas(1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(2,DER_FUNC2D,i,j) =  0.5_DP*dx-0.5_DP*dy
              Dbas(3,DER_FUNC2D,i,j) =  0.5_DP*dx+0.5_DP*dy
              Dbas(4,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(5,DER_FUNC2D,i,j) = 1.0_DP - dx

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas(1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(3,DER_FUNC2D,i,j) =  0.5_DP*dx+0.5_DP*dy 
              Dbas(4,DER_FUNC2D,i,j) = -0.5_DP*dx+0.5_DP*dy 
              Dbas(5,DER_FUNC2D,i,j) = 1.0_DP - dy

            else
              ! We are in T4

              Dbas(1,DER_FUNC2D,i,j) = -0.5_DP*dx-0.5_DP*dy
              Dbas(2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas(4,DER_FUNC2D,i,j) = -0.5_DP*dx+0.5_DP*dy
              Dbas(5,DER_FUNC2D,i,j) = 1.0_DP + dx

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,ddet,dx,dy,a11,a12,a21,a22)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)
          a11 = reval%p_Djac(1,i,j)
          a21 = reval%p_Djac(2,i,j)
          a12 = reval%p_Djac(3,i,j)
          a22 = reval%p_Djac(4,i,j)

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas(1,DER_DERIV2D_X,i,j) = ((-0.5_DP)*a22 - (-0.5_DP)*a21)*ddet
              Dbas(2,DER_DERIV2D_X,i,j) = (( 0.5_DP)*a22 - (-0.5_DP)*a21)*ddet
              Dbas(3,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(4,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(5,DER_DERIV2D_X,i,j) = (-(1.0_DP)*a21)*ddet

              Dbas(1,DER_DERIV2D_Y,i,j) = (-(-0.5_DP)*a12 + (-0.5_DP)*a11)*ddet
              Dbas(2,DER_DERIV2D_Y,i,j) = (-( 0.5_DP)*a12 + (-0.5_DP)*a11)*ddet
              Dbas(3,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(4,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(5,DER_DERIV2D_Y,i,j) = ((1.0_DP)*a11)*ddet

            else
              ! We are in T2

              Dbas(1,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(2,DER_DERIV2D_X,i,j) = (( 0.5_DP)*a22 - (-0.5_DP)*a21)*ddet
              Dbas(3,DER_DERIV2D_X,i,j) = (( 0.5_DP)*a22 - ( 0.5_DP)*a21)*ddet
              Dbas(4,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(5,DER_DERIV2D_X,i,j) = ((-1.0_DP)*a22)*ddet

              Dbas(1,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(2,DER_DERIV2D_Y,i,j) = (-( 0.5_DP)*a12 + (-0.5_DP)*a11)*ddet
              Dbas(3,DER_DERIV2D_Y,i,j) = (-( 0.5_DP)*a12 + ( 0.5_DP)*a11)*ddet
              Dbas(4,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(5,DER_DERIV2D_Y,i,j) = (-(-1.0_DP)*a12)*ddet

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas(1,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(2,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(3,DER_DERIV2D_X,i,j) = (( 0.5_DP)*a22 - (0.5_DP)*a21)*ddet
              Dbas(4,DER_DERIV2D_X,i,j) = ((-0.5_DP)*a22 - (0.5_DP)*a21)*ddet
              Dbas(5,DER_DERIV2D_X,i,j) = (-(-1.0_DP)*a21)*ddet

              Dbas(1,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(2,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(3,DER_DERIV2D_Y,i,j) = (-( 0.5_DP)*a12 + (0.5_DP)*a11)*ddet
              Dbas(4,DER_DERIV2D_Y,i,j) = (-(-0.5_DP)*a12 + (0.5_DP)*a11)*ddet
              Dbas(5,DER_DERIV2D_Y,i,j) = ((-1.0_DP)*a11)*ddet

            else
              ! We are in T4

              Dbas(1,DER_DERIV2D_X,i,j) = ((-0.5_DP)*a22 - (-0.5_DP)*a21)*ddet
              Dbas(2,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(3,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(4,DER_DERIV2D_X,i,j) = ((-0.5_DP)*a22 - (0.5_DP)*a21)*ddet
              Dbas(5,DER_DERIV2D_X,i,j) = ((1.0_DP)*a22)*ddet

              Dbas(1,DER_DERIV2D_Y,i,j) = (-(-0.5_DP)*a12 + (-0.5_DP)*a11)*ddet
              Dbas(2,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(3,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(4,DER_DERIV2D_Y,i,j) = (-(-0.5_DP)*a12 + (0.5_DP)*a11)*ddet
              Dbas(5,DER_DERIV2D_Y,i,j) = (-(1.0_DP)*a12)*ddet

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_QPW4P1T_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! This element is piecewisely defined by 4 triangles in one
  ! quad. There are 8 DOFs: the four edge midpoints plus the four
  ! midpoints along the intersected diagonals:
  !
  ! -1,1              1,1       a3                a2
  !   +-------------+             +------*------+
  !   | \    T3   / |             | \         / |
  !   |   \     /   |             |   *     *   |
  !   |     \ /     |             |     \ /     |
  !   | T4   X   T2 |     ->      *      X MP   *
  !   |     / \     |             |     / \     |
  !   |   /     \   |             |   *     *   |
  !   | /    T1   \ |             | /         \ |
  !   +-------------+             +------*------+
  ! -1,-1             1,-1      a0                a1
  !
  ! Therefore, for every point we have to decide in which
  ! sub-triangle we are.
  !
  ! The point ia in one of the subtriangles T1, T2, T3 or T4.
  ! Calculate the line a0->a2 and chech whether it is left
  ! or right. Calculate the line a1->a3 and check whether it is
  ! left or right. That way, we know in which subtriangle the
  ! point is.

  ! Local variables
  real(DP) :: ddet,dx,dy,a11,a12,a21,a22
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas(1,DER_FUNC2D,i,j) = -1.0_DP-2.0_DP*dy
              Dbas(2,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(3,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(4,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(5,DER_FUNC2D,i,j) =  1.0_DP+dx+dy
              Dbas(6,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(7,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(8,DER_FUNC2D,i,j) =  1.0_DP-dx+dy

            else
              ! We are in T2

              Dbas(1,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(2,DER_FUNC2D,i,j) = -1.0_DP+2.0_DP*dx
              Dbas(3,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(4,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(5,DER_FUNC2D,i,j) =  1.0_DP-dx-dy
              Dbas(6,DER_FUNC2D,i,j) =  1.0_DP-dx+dy
              Dbas(7,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(8,DER_FUNC2D,i,j) =  0.0_DP

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas(1,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(2,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(3,DER_FUNC2D,i,j) = -1.0_DP+2.0_DP*dy
              Dbas(4,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(5,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(6,DER_FUNC2D,i,j) =  1.0_DP+dx-dy
              Dbas(7,DER_FUNC2D,i,j) =  1.0_DP-dx-dy
              Dbas(8,DER_FUNC2D,i,j) =  0.0_DP

            else
              ! We are in T4

              Dbas(1,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(2,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(3,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(4,DER_FUNC2D,i,j) = -1.0_DP-2.0_DP*dx
              Dbas(5,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(6,DER_FUNC2D,i,j) =  0.0_DP
              Dbas(7,DER_FUNC2D,i,j) =  1.0_DP+dx+dy
              Dbas(8,DER_FUNC2D,i,j) =  1.0_DP+dx-dy

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,ddet,dx,dy,a11,a12,a21,a22)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)
          a11 = reval%p_Djac(1,i,j)
          a21 = reval%p_Djac(2,i,j)
          a12 = reval%p_Djac(3,i,j)
          a22 = reval%p_Djac(4,i,j)

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas(1,DER_DERIV2D_X,i,j) = -(-2.0_DP)*a21*ddet
              Dbas(2,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(3,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(4,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(5,DER_DERIV2D_X,i,j) =  (a22-a21)*ddet
              Dbas(6,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(7,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(8,DER_DERIV2D_X,i,j) =  (-a22-a21)*ddet

              Dbas(1,DER_DERIV2D_Y,i,j) = -2.0_DP*a11*ddet
              Dbas(2,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(3,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(4,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(5,DER_DERIV2D_Y,i,j) =  (-a12+a11)*ddet
              Dbas(6,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(7,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(8,DER_DERIV2D_Y,i,j) =  (a12+a11)*ddet

            else
              ! We are in T2

              Dbas(1,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(2,DER_DERIV2D_X,i,j) =  2.0_DP*a22*ddet
              Dbas(3,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(4,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(5,DER_DERIV2D_X,i,j) =  (-a22+a21)*ddet 
              Dbas(6,DER_DERIV2D_X,i,j) =  (-a22-a21)*ddet 
              Dbas(7,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(8,DER_DERIV2D_X,i,j) =  0.0_DP

              Dbas(1,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(2,DER_DERIV2D_Y,i,j) = -2.0_DP*a12*ddet
              Dbas(3,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(4,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(5,DER_DERIV2D_Y,i,j) =  (a12-a11)*ddet
              Dbas(6,DER_DERIV2D_Y,i,j) =  (a12+a11)*ddet
              Dbas(7,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(8,DER_DERIV2D_Y,i,j) =  0.0_DP

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas(1,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(2,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(3,DER_DERIV2D_X,i,j) = -2.0_DP*a21*ddet
              Dbas(4,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(5,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(6,DER_DERIV2D_X,i,j) =  (a22+a21)*ddet
              Dbas(7,DER_DERIV2D_X,i,j) =  (-a22+a21)*ddet
              Dbas(8,DER_DERIV2D_X,i,j) =  0.0_DP

              Dbas(1,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(2,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(3,DER_DERIV2D_Y,i,j) =  2.0_DP*a11*ddet
              Dbas(4,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(5,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(6,DER_DERIV2D_Y,i,j) =  (-a12-a11)*ddet
              Dbas(7,DER_DERIV2D_Y,i,j) =  (a12-a11)*ddet
              Dbas(8,DER_DERIV2D_Y,i,j) =  0.0_DP

            else
              ! We are in T4

              Dbas(1,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(2,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(3,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(4,DER_DERIV2D_X,i,j) = -2.0_DP*a22*ddet
              Dbas(5,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(6,DER_DERIV2D_X,i,j) =  0.0_DP
              Dbas(7,DER_DERIV2D_X,i,j) =  (a22-a21)*ddet
              Dbas(8,DER_DERIV2D_X,i,j) =  (a22+a21)*ddet

              Dbas(1,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(2,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(3,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(4,DER_DERIV2D_Y,i,j) =  2.0_DP*a12*ddet
              Dbas(5,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(6,DER_DERIV2D_Y,i,j) =  0.0_DP
              Dbas(7,DER_DERIV2D_Y,i,j) =  (-a12+a11)*ddet
              Dbas(8,DER_DERIV2D_Y,i,j) =  (-a12-a11)*ddet

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine
  
  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_QPW4DCP1_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! This element is piecewisely defined by 4 triangles in one
  ! quad. There are 5 DOFs: the four corners plus the element
  ! midpoint:
  !
  ! -1,1              1,1       a3                a2
  !   +-------------+             +-------------+
  !   | \    T3   / |             | \         / |
  !   |   \     /   |             |   \     /   |
  !   |     \ /     |             |     \ /     |
  !   | T4   X   T2 |     ->      |      X MP   |
  !   |     / \     |             |     / \     |
  !   |   /     \   |             |   /     \   |
  !   | /    T1   \ |             | /         \ |
  !   +-------------+             +-------------+
  ! -1,-1             1,-1      a0                a1
  !
  ! Therefore, for every point we have to decide in which
  ! sub-triangle we are.
  !
  ! The point ia in one of the subtriangles T1, T2, T3 or T4.
  ! Calculate the line a0->a2 and chech whether it is left
  ! or right. Calculate the line a1->a3 and check whether it is
  ! left or right. That way, we know in which subtriangle the
  ! point is.

  ! Local variables
  real(DP) :: ddet,dx,dy
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.
              Dbas( 1,DER_FUNC2D,i,j) = -0.5_DP*dx-0.5_DP*dy
              Dbas( 2,DER_FUNC2D,i,j) =  0.5_DP*dx-0.5_DP*dy
              Dbas( 3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 4,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 5,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 6,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 7,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 8,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 9,DER_FUNC2D,i,j) = 1.0_DP + dy
              Dbas(10,DER_FUNC2D,i,j) = 1.0_DP + dy
              Dbas(11,DER_FUNC2D,i,j) = 1.0_DP + dy

            else
              ! We are in T2
              Dbas( 1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 3,DER_FUNC2D,i,j) =  0.5_DP*dx-0.5_DP*dy
              Dbas( 4,DER_FUNC2D,i,j) =  0.5_DP*dx+0.5_DP*dy
              Dbas( 5,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 6,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 7,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 8,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 9,DER_FUNC2D,i,j) = 1.0_DP - dx
              Dbas(10,DER_FUNC2D,i,j) = 1.0_DP - dx
              Dbas(11,DER_FUNC2D,i,j) = -1.0_DP + dx

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3
              Dbas( 1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 4,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 5,DER_FUNC2D,i,j) =  0.5_DP*dx+0.5_DP*dy 
              Dbas( 6,DER_FUNC2D,i,j) = -0.5_DP*dx+0.5_DP*dy 
              Dbas( 7,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 8,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 9,DER_FUNC2D,i,j) = 1.0_DP - dy
              Dbas(10,DER_FUNC2D,i,j) = -1.0_DP + dy
              Dbas(11,DER_FUNC2D,i,j) = -1.0_DP + dy

            else
              ! We are in T4
              Dbas( 1,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 2,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 3,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 4,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 5,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 6,DER_FUNC2D,i,j) = 0.0_DP
              Dbas( 7,DER_FUNC2D,i,j) = -0.5_DP*dx+0.5_DP*dy
              Dbas( 8,DER_FUNC2D,i,j) = -0.5_DP*dx-0.5_DP*dy
              Dbas( 9,DER_FUNC2D,i,j) = 1.0_DP + dx
              Dbas(10,DER_FUNC2D,i,j) = -1.0_DP - dx
              Dbas(11,DER_FUNC2D,i,j) = 1.0_DP + dx

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas( 1,DER_DERIV2D_X,i,j) = 0.5_DP*(-reval%p_Djac(4,i,j) + reval%p_Djac(2,i,j))*ddet
              Dbas( 2,DER_DERIV2D_X,i,j) = 0.5_DP*( reval%p_Djac(4,i,j) + reval%p_Djac(2,i,j))*ddet
              Dbas( 3,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_X,i,j) = -reval%p_Djac(2,i,j)*ddet
              Dbas(10,DER_DERIV2D_X,i,j) = -reval%p_Djac(2,i,j)*ddet
              Dbas(11,DER_DERIV2D_X,i,j) = -reval%p_Djac(2,i,j)*ddet

              Dbas( 1,DER_DERIV2D_Y,i,j) = 0.5_DP*( reval%p_Djac(3,i,j) - reval%p_Djac(1,i,j))*ddet
              Dbas( 2,DER_DERIV2D_Y,i,j) = 0.5_DP*(-reval%p_Djac(3,i,j) - reval%p_Djac(1,i,j))*ddet
              Dbas( 3,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_Y,i,j) = reval%p_Djac(1,i,j)*ddet
              Dbas(10,DER_DERIV2D_Y,i,j) = reval%p_Djac(1,i,j)*ddet
              Dbas(11,DER_DERIV2D_Y,i,j) = reval%p_Djac(1,i,j)*ddet

            else
              ! We are in T2

              Dbas( 1,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_X,i,j) = 0.5_DP*(reval%p_Djac(4,i,j) + reval%p_Djac(2,i,j))*ddet
              Dbas( 4,DER_DERIV2D_X,i,j) = 0.5_DP*(reval%p_Djac(4,i,j) - reval%p_Djac(2,i,j))*ddet
              Dbas( 5,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_X,i,j) = -reval%p_Djac(4,i,j)*ddet
              Dbas(10,DER_DERIV2D_X,i,j) = -reval%p_Djac(4,i,j)*ddet
              Dbas(11,DER_DERIV2D_X,i,j) =  reval%p_Djac(4,i,j)*ddet

              Dbas( 1,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_Y,i,j) = 0.5_DP*(-reval%p_Djac(3,i,j) - reval%p_Djac(1,i,j))*ddet
              Dbas( 4,DER_DERIV2D_Y,i,j) = 0.5_DP*(-reval%p_Djac(3,i,j) + reval%p_Djac(1,i,j))*ddet
              Dbas( 5,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_Y,i,j) =  reval%p_Djac(3,i,j)*ddet
              Dbas(10,DER_DERIV2D_Y,i,j) =  reval%p_Djac(3,i,j)*ddet
              Dbas(11,DER_DERIV2D_Y,i,j) = -reval%p_Djac(3,i,j)*ddet

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas( 1,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_X,i,j) = 0.5_DP*( reval%p_Djac(4,i,j) - reval%p_Djac(2,i,j))*ddet
              Dbas( 6,DER_DERIV2D_X,i,j) = 0.5_DP*(-reval%p_Djac(4,i,j) - reval%p_Djac(2,i,j))*ddet
              Dbas( 7,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_X,i,j) =  reval%p_Djac(2,i,j)*ddet
              Dbas(10,DER_DERIV2D_X,i,j) = -reval%p_Djac(2,i,j)*ddet
              Dbas(11,DER_DERIV2D_X,i,j) = -reval%p_Djac(2,i,j)*ddet

              Dbas( 1,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_Y,i,j) = 0.5_DP*(-reval%p_Djac(3,i,j) + reval%p_Djac(1,i,j))*ddet
              Dbas( 6,DER_DERIV2D_Y,i,j) = 0.5_DP*( reval%p_Djac(3,i,j) + reval%p_Djac(1,i,j))*ddet
              Dbas( 7,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_Y,i,j) = -reval%p_Djac(1,i,j)*ddet
              Dbas(10,DER_DERIV2D_Y,i,j) =  reval%p_Djac(1,i,j)*ddet
              Dbas(11,DER_DERIV2D_Y,i,j) =  reval%p_Djac(1,i,j)*ddet

            else
              ! We are in T4
              Dbas( 1,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_X,i,j) = 0.5_DP*(-reval%p_Djac(4,i,j) - reval%p_Djac(2,i,j))*ddet
              Dbas( 8,DER_DERIV2D_X,i,j) = 0.5_DP*(-reval%p_Djac(4,i,j) + reval%p_Djac(2,i,j))*ddet
              Dbas( 9,DER_DERIV2D_X,i,j) =  reval%p_Djac(4,i,j)*ddet
              Dbas(10,DER_DERIV2D_X,i,j) = -reval%p_Djac(4,i,j)*ddet
              Dbas(11,DER_DERIV2D_X,i,j) =  reval%p_Djac(4,i,j)*ddet

              Dbas( 1,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_Y,i,j) = 0.5_DP*(reval%p_Djac(3,i,j) + reval%p_Djac(1,i,j))*ddet
              Dbas( 8,DER_DERIV2D_Y,i,j) = 0.5_DP*(reval%p_Djac(3,i,j) - reval%p_Djac(1,i,j))*ddet
              Dbas( 9,DER_DERIV2D_Y,i,j) = -reval%p_Djac(3,i,j)*ddet
              Dbas(10,DER_DERIV2D_Y,i,j) =  reval%p_Djac(3,i,j)*ddet
              Dbas(11,DER_DERIV2D_Y,i,j) = -reval%p_Djac(3,i,j)*ddet

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_QPW4P2_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! This element is piecewisely defined by 4 triangles in one
  ! quad. There are 13 local DOFs: 
  !  - the four corners
  !  - the four edge midpoints
  !  - the element center
  !  - the four edge midpoints of the inner edges
  ! On every triangle, we have a P2 space.
  !
  !  -1,1                           1,1
  !    (4)-----------(7)-----------(3)
  !     |`.                       ,'|
  !     |  `.                   ,'  |
  !     |    `.       T3      ,'    |
  !     |     (13)         (12)     |
  !     |        `.       ,'        |
  !     |          `.   ,'          |
  !    (8)   T4      (9)      T2   (6)
  !     |          ,'   `.          |
  !     |        ,'       `.        |
  !     |     (10)         (11)     |
  !     |    ,'       T1      `.    |
  !     |  ,'                   `.  |
  !     |,'                       `.|
  !    (1)-----------(5)-----------(2)
  ! -1,-1                           1,-1  
  !
  ! For every point we have to decide in which
  ! sub-triangle we are.
  !
  ! The point ia in one of the subtriangles T1, T2, T3 or T4.
  ! Calculate the line a0->a2 and chech whether it is left
  ! or right. Calculate the line a1->a3 and check whether it is
  ! left or right. That way, we know in which subtriangle the
  ! point is.

  ! Local variables
  real(DP) :: ddet,dx,dy,a11,a12,a21,a22
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas( 1,DER_FUNC2D,i,j) = dx/2+dy/2+dx**2/2+dx*dy+dy**2/2
              Dbas( 2,DER_FUNC2D,i,j) = -dx/2+dy/2+dx**2/2-dx*dy+dy**2/2
              Dbas( 3,DER_FUNC2D,i,j) = 0
              Dbas( 4,DER_FUNC2D,i,j) = 0
              Dbas( 5,DER_FUNC2D,i,j) = -dx**2+dy**2
              Dbas( 6,DER_FUNC2D,i,j) = 0
              Dbas( 7,DER_FUNC2D,i,j) = 0
              Dbas( 8,DER_FUNC2D,i,j) = 0
              Dbas( 9,DER_FUNC2D,i,j) = 1+3*dy+2*dy**2
              Dbas(10,DER_FUNC2D,i,j) = -2*dx-2*dy-2*dx*dy-2*dy**2
              Dbas(11,DER_FUNC2D,i,j) = 2*dx-2*dy+2*dx*dy-2*dy**2
              Dbas(12,DER_FUNC2D,i,j) = 0
              Dbas(13,DER_FUNC2D,i,j) = 0

            else
              ! We are in T2

              Dbas( 1,DER_FUNC2D,i,j) = 0
              Dbas( 2,DER_FUNC2D,i,j) = -dx/2+dy/2+dx**2/2-dx*dy+dy**2/2
              Dbas( 3,DER_FUNC2D,i,j) = -dx/2-dy/2+dx**2/2+dx*dy+dy**2/2
              Dbas( 4,DER_FUNC2D,i,j) = 0
              Dbas( 5,DER_FUNC2D,i,j) = 0
              Dbas( 6,DER_FUNC2D,i,j) = dx**2-dy**2
              Dbas( 7,DER_FUNC2D,i,j) = 0
              Dbas( 8,DER_FUNC2D,i,j) = 0
              Dbas( 9,DER_FUNC2D,i,j) = 1-3*dx+2*dx**2
              Dbas(10,DER_FUNC2D,i,j) = 0
              Dbas(11,DER_FUNC2D,i,j) = 2*dx-2*dy-2*dx**2+2*dx*dy
              Dbas(12,DER_FUNC2D,i,j) = 2*dx+2*dy-2*dx**2-2*dx*dy
              Dbas(13,DER_FUNC2D,i,j) = 0

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas( 1,DER_FUNC2D,i,j) = 0
              Dbas( 2,DER_FUNC2D,i,j) = 0
              Dbas( 3,DER_FUNC2D,i,j) = -dx/2-dy/2+dx**2/2+dx*dy+dy**2/2
              Dbas( 4,DER_FUNC2D,i,j) = dx/2-dy/2+dx**2/2-dx*dy+dy**2/2
              Dbas( 5,DER_FUNC2D,i,j) = 0
              Dbas( 6,DER_FUNC2D,i,j) = 0
              Dbas( 7,DER_FUNC2D,i,j) = -dx**2+dy**2
              Dbas( 8,DER_FUNC2D,i,j) = 0
              Dbas( 9,DER_FUNC2D,i,j) = 1-3*dy+2*dy**2
              Dbas(10,DER_FUNC2D,i,j) = 0
              Dbas(11,DER_FUNC2D,i,j) = 0
              Dbas(12,DER_FUNC2D,i,j) = 2*dx+2*dy-2*dx*dy-2*dy**2
              Dbas(13,DER_FUNC2D,i,j) = -2*dx+2*dy+2*dx*dy-2*dy**2

            else
              ! We are in T4

              Dbas( 1,DER_FUNC2D,i,j) = dx/2+dy/2+dx**2/2+dx*dy+dy**2/2
              Dbas( 2,DER_FUNC2D,i,j) = 0
              Dbas( 3,DER_FUNC2D,i,j) = 0
              Dbas( 4,DER_FUNC2D,i,j) = dx/2-dy/2+dx**2/2-dx*dy+dy**2/2
              Dbas( 5,DER_FUNC2D,i,j) = 0
              Dbas( 6,DER_FUNC2D,i,j) = 0
              Dbas( 7,DER_FUNC2D,i,j) = 0
              Dbas( 8,DER_FUNC2D,i,j) = dx**2-dy**2
              Dbas( 9,DER_FUNC2D,i,j) = 1+3*dx+2*dx**2
              Dbas(10,DER_FUNC2D,i,j) = -2*dx-2*dy-2*dx**2-2*dx*dy
              Dbas(11,DER_FUNC2D,i,j) = 0
              Dbas(12,DER_FUNC2D,i,j) = 0
              Dbas(13,DER_FUNC2D,i,j) = -2*dx+2*dy-2*dx**2+2*dx*dy

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,ddet,dx,dy,a11,a12,a21,a22)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          a11 = reval%p_Djac(1,i,j)
          a21 = reval%p_Djac(2,i,j)
          a12 = reval%p_Djac(3,i,j)
          a22 = reval%p_Djac(4,i,j)

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! Figure out in which triangle we are
          if (dy .le. dx) then
            ! We are in T1 or T2.
            if (dy .le. -dx) then
              ! We are in T1.

              Dbas( 1,DER_DERIV2D_X,i,j) = ddet*(a22-a21)*(2*dy+1+2*dx)/2
              Dbas( 2,DER_DERIV2D_X,i,j) = ddet*(a22+a21)*(-2*dy-1+2*dx)/2
              Dbas( 3,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_X,i,j) = -2*ddet*(a22*dx+a21*dy)
              Dbas( 6,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_X,i,j) = -a21*ddet*(3+4*dy)
              Dbas(10,DER_DERIV2D_X,i,j) = 2*ddet*(-a22-a22*dy+a21+a21*dx+2*a21*dy)
              Dbas(11,DER_DERIV2D_X,i,j) = -2*ddet*(-a22-dy*a22-a21+a21*dx-2*a21*dy)
              Dbas(12,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(13,DER_DERIV2D_X,i,j) = 0.0_DP

              Dbas( 1,DER_DERIV2D_Y,i,j) = ddet*(-a12+a11)*(1+2*dx+2*dy)/2
              Dbas( 2,DER_DERIV2D_Y,i,j) = -ddet*(a12+a11)*(-1+2*dx-2*dy)/2
              Dbas( 3,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_Y,i,j) = 2*ddet*(a12*dx+a11*dy)
              Dbas( 6,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_Y,i,j) = a11*ddet*(3+4*dy)
              Dbas(10,DER_DERIV2D_Y,i,j) = -2*ddet*(-a12-dy*a12+a11+a11*dx+2*a11*dy)
              Dbas(11,DER_DERIV2D_Y,i,j) = 2*ddet*(-a12-dy*a12-a11+a11*dx-2*a11*dy)
              Dbas(12,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(13,DER_DERIV2D_Y,i,j) = 0.0_DP

            else
              ! We are in T2

              Dbas( 1,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_X,i,j) = ddet*(a22+a21)*(-2*dy-1+2*dx)/2
              Dbas( 3,DER_DERIV2D_X,i,j) = ddet*(a22-a21)*(2*dy-1+2*dx)/2
              Dbas( 4,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_X,i,j) = 2*ddet*(a22*dx+dy*a21)
              Dbas( 7,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_X,i,j) = a22*ddet*(-3+4*dx)
              Dbas(10,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(11,DER_DERIV2D_X,i,j) = -2*ddet*(-a22+2*a22*dx-a22*dy-a21+a21*dx)
              Dbas(12,DER_DERIV2D_X,i,j) = -2*ddet*(-a22+2*a22*dx+a22*dy+a21-a21*dx)
              Dbas(13,DER_DERIV2D_X,i,j) = 0.0_DP

              Dbas( 1,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_Y,i,j) = -ddet*(a12+a11)*(-1+2*dx-2*dy)/2
              Dbas( 3,DER_DERIV2D_Y,i,j) = ddet*(-a12+a11)*(-1+2*dx+2*dy)/2
              Dbas( 4,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 5,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_Y,i,j) = -2*ddet*(a12*dx+a11*dy)
              Dbas( 7,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_Y,i,j) = -a12*ddet*(-3+4*dx)
              Dbas(10,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(11,DER_DERIV2D_Y,i,j) = 2*ddet*(-a12+2*a12*dx-a12*dy-a11+a11*dx)
              Dbas(12,DER_DERIV2D_Y,i,j) = -2*ddet*(a12-2*a12*dx-a12*dy-a11+a11*dx)
              Dbas(13,DER_DERIV2D_Y,i,j) = 0.0_DP

            end if
          else
            ! We are in T3 or T4
            if (dy .gt. -dx) then
              ! We are in T3

              Dbas( 1,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_X,i,j) = ddet*(a22-a21)*(2*dy-1+2*dx)/2
              Dbas( 4,DER_DERIV2D_X,i,j) = ddet*(a22+a21)*(-2*dy+1+2*dx)/2
              Dbas( 5,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_X,i,j) = -2*ddet*(a22*dx+dy*a21)
              Dbas( 8,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_X,i,j) = -a21*ddet*(-3+4*dy)
              Dbas(10,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(11,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(12,DER_DERIV2D_X,i,j) = 2*ddet*(a22-dy*a22-a21+a21*dx+2*dy*a21)
              Dbas(13,DER_DERIV2D_X,i,j) = -2*ddet*(a22-dy*a22+a21+a21*dx-2*dy*a21)

              Dbas( 1,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 2,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_Y,i,j) = ddet*(-a12+a11)*(-1+2*dx+2*dy)/2
              Dbas( 4,DER_DERIV2D_Y,i,j) = -ddet*(a12+a11)*(1+2*dx-2*dy)/2
              Dbas( 5,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_Y,i,j) = 2*ddet*(a12*dx+a11*dy)
              Dbas( 8,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 9,DER_DERIV2D_Y,i,j) = a11*ddet*(-3+4*dy)
              Dbas(10,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(11,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(12,DER_DERIV2D_Y,i,j) = -2*ddet*(a12-a12*dy-a11+a11*dx+2*a11*dy)
              Dbas(13,DER_DERIV2D_Y,i,j) = 2*ddet*(a12-dy*a12+a11+a11*dx-2*a11*dy)

            else
              ! We are in T4

              Dbas( 1,DER_DERIV2D_X,i,j) = ddet*(a22-a21)*(2*dy+1+2*dx)/2
              Dbas( 2,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_X,i,j) = ddet*(a22+a21)*(-2*dy+1+2*dx)/2
              Dbas( 5,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_X,i,j) = 2*ddet*(a22*dx+a21*dy)
              Dbas( 9,DER_DERIV2D_X,i,j) = a22*ddet*(3+4*dx)
              Dbas(10,DER_DERIV2D_X,i,j) = -2*ddet*(a22+2*a22*dx+a22*dy-a21-a21*dx)
              Dbas(11,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(12,DER_DERIV2D_X,i,j) = 0.0_DP
              Dbas(13,DER_DERIV2D_X,i,j) = -2*ddet*(a22+2*a22*dx-a22*dy+a21+a21*dx)

              Dbas( 1,DER_DERIV2D_Y,i,j) = ddet*(-a12+a11)*(1+2*dx+2*dy)/2
              Dbas( 2,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 3,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 4,DER_DERIV2D_Y,i,j) = -ddet*(a12+a11)*(1+2*dx-2*dy)/2
              Dbas( 5,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 6,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 7,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas( 8,DER_DERIV2D_Y,i,j) = -2*ddet*(a12*dx+dy*a11)
              Dbas( 9,DER_DERIV2D_Y,i,j) = -a12*ddet*(3+4*dx)
              Dbas(10,DER_DERIV2D_Y,i,j) = -2*ddet*(-a12-2*a12*dx-a12*dy+a11+a11*dx)
              Dbas(11,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(12,DER_DERIV2D_Y,i,j) = 0.0_DP
              Dbas(13,DER_DERIV2D_Y,i,j) = 2*ddet*(a12+2*a12*dx-a12*dy+a11+a11*dx)

            end if
          end if

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_E030_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The E030_2D element is specified by four polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x^2 - y^2 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,4:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))) d(t) = kronecker(i,j) * |ej|
  !     <==>
  !     Int_ej (Pi(t)) d(t) = kronecker(i,j) * |ej|
  !   }
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  !  P1(x,y) = -3/8 * (x^2 - y^2)  -  1/2 * y  +  1/4
  !  P2(x,y) =  3/8 * (x^2 - y^2)  +  1/2 * x  +  1/4
  !  P3(x,y) = -3/8 * (x^2 - y^2)  +  1/2 * y  +  1/4
  !  P4(x,y) =  3/8 * (x^2 - y^2)  -  1/2 * x  +  1/4

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 4

  ! Local variables
  real(DP) :: ddet,dx,dy,dt
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy,dt)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)
          dt = 0.375_DP * (dx*dx - dy*dy)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = -dt - 0.5_DP*dy + 0.25_DP
          Dbas(2,DER_FUNC2D,i,j) =  dt + 0.5_DP*dx + 0.25_DP
          Dbas(3,DER_FUNC2D,i,j) = -dt + 0.5_DP*dy + 0.25_DP
          Dbas(4,DER_FUNC2D,i,j) =  dt - 0.5_DP*dx + 0.25_DP

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) = -0.75_DP*dx
          DrefDer(2,1) =  0.75_DP*dx + 0.5_DP
          DrefDer(3,1) = -0.75_DP*dx
          DrefDer(4,1) =  0.75_DP*dx - 0.5_DP
          ! Y-derivatives
          DrefDer(1,2) = -0.75_DP*dy - 0.5_DP
          DrefDer(2,2) =  0.75_DP*dy
          DrefDer(3,2) = -0.75_DP*dy + 0.5_DP
          DrefDer(4,2) =  0.75_DP*dy

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_EB30_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The EB30_2D element is specified by five polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2 - y^2 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,5:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))) d(t) = kronecker(i,j) * |ej|
  !     <==>
  !     Int_ej (Pi(t)) d(t) = kronecker(i,j) * |ej|
  !   }
  !   Int_T (p_i(x,y)*L1(x)*L1(y)) d(x,y) = kronecker(i,5) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  ! T being the quadrilateral
  ! |T| being the area of the quadrilateral
  ! L1 being the first Legendre-Polynomial:
  ! L1(x) := x
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  !  P1(x,y) = -3/8 * (x^2 - y^2)  -  1/2 * y  +  1/4
  !  P2(x,y) =  3/8 * (x^2 - y^2)  +  1/2 * x  +  1/4
  !  P3(x,y) = -3/8 * (x^2 - y^2)  +  1/2 * y  +  1/4
  !  P4(x,y) =  3/8 * (x^2 - y^2)  -  1/2 * x  +  1/4
  !  P5(x,y) =  9*x*y

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 5

  ! Local variables
  real(DP) :: ddet,dx,dy,dt
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy,dt)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)
          dt = 0.375_DP * (dx*dx - dy*dy)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = -dt - 0.5_DP*dy + 0.25_DP
          Dbas(2,DER_FUNC2D,i,j) =  dt + 0.5_DP*dx + 0.25_DP
          Dbas(3,DER_FUNC2D,i,j) = -dt + 0.5_DP*dy + 0.25_DP
          Dbas(4,DER_FUNC2D,i,j) =  dt - 0.5_DP*dx + 0.25_DP
          Dbas(5,DER_FUNC2D,i,j) = 9.0_DP*dx*dy

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) = -0.75_DP*dx
          DrefDer(2,1) =  0.75_DP*dx + 0.5_DP
          DrefDer(3,1) = -0.75_DP*dx
          DrefDer(4,1) =  0.75_DP*dx - 0.5_DP
          DrefDer(5,1) =  9.0_DP*dy
          ! Y-derivatives
          DrefDer(1,2) = -0.75_DP*dy - 0.5_DP
          DrefDer(2,2) =  0.75_DP*dy
          DrefDer(3,2) = -0.75_DP*dy + 0.5_DP
          DrefDer(4,2) =  0.75_DP*dy
          DrefDer(5,2) =  9.0_DP*dx

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_Q1TBNP_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The Q1TBNP_2D element is specified by five polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2 - y^2 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,5:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))) d(t) = kronecker(i,j) * |ej|
  !     <==>
  !     Int_ej (Pi(t)) d(t) = kronecker(i,j) * |ej|
  !   }
  !   Int_T (Pi(x,y)*L1(x)*L1(y)) d(x,y) = kronecker(i,5) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t



  ! Parameter: Number of local basis functions
  integer, parameter :: NBAS = 5

  ! Parameter: Number of cubature points for 1D edge integration
  integer, parameter :: NCUB1D = 2
  !integer, parameter :: NCUB1D = 5

  ! Parameter: Number of cubature points for 2D quad integration
  integer, parameter :: NCUB2D = NCUB1D**2

  ! 1D edge cubature rule point coordinates and weights
  real(DP), dimension(NCUB1D) :: DcubPts1D
  real(DP), dimension(NCUB1D) :: DcubOmega1D

  ! 2D quad cubature rule point coordinates and weights
  real(DP), dimension(NDIM2D, NCUB2D) :: DcubPts2D
  real(DP), dimension(NCUB2D) :: DcubOmega2D

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert

  ! Local mapped 1D cubature point coordinates and integration weights
  real(DP), dimension(NDIM2D, NCUB1D, 4) :: DedgePoints
  real(DP), dimension(NCUB1D, 4) :: DedgeWeights
  real(DP), dimension(4) :: DedgeLen

  ! Local mapped 2D cubature point coordinates and jacobian determinants
  real(DP), dimension(NDIM2D, NCUB2D) :: DquadPoints
  real(DP), dimension(NCUB2D) :: DquadWeights
  real(DP) :: dquadArea

  ! temporary variables for trafo call (will be removed later)
  real(DP), dimension(8) :: DjacPrep
  real(DP), dimension(4) :: DjacTrafo

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr
  real(DP) :: ddets

  ! other local variables
  integer :: i,j,k,iel, ipt
  real(DP), dimension(NBAS,NBAS) :: Da, Dc
  real(DP) :: dx,dy,dt,derx,dery
  logical :: bsuccess


    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Step 0: Set up 1D and 2D cubature rules
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Set up a 2-point Gauss rule for 1D
    DcubPts1D(1) = -sqrt(1.0_DP / 3.0_DP)
    DcubPts1D(2) =  sqrt(1.0_DP / 3.0_DP)
    DcubOmega1D(1) = 1.0_DP
    DcubOmega1D(2) = 1.0_DP

!    ! !!! DEBUG: 5-point Gauss rule !!!
!    dt = 2.0_DP*sqrt(10.0_DP / 7.0_DP)
!    DcubPts1D(1) = -sqrt(5.0_DP + dt) / 3.0_DP
!    DcubPts1D(2) = -sqrt(5.0_DP - dt) / 3.0_DP
!    DcubPts1D(3) = 0.0_DP
!    DcubPts1D(4) =  sqrt(5.0_DP - dt) / 3.0_DP
!    DcubPts1D(5) =  sqrt(5.0_DP + dt) / 3.0_DP
!    dt = 13.0_DP*sqrt(70.0_DP)
!    DcubOmega1D(1) = (322.0_DP - dt) / 900.0_DP
!    DcubOmega1D(2) = (322.0_DP + dt) / 900.0_DP
!    DcubOmega1D(3) = 128.0_DP / 225.0_DP
!    DcubOmega1D(4) = (322.0_DP + dt) / 900.0_DP
!    DcubOmega1D(5) = (322.0_DP - dt) / 900.0_DP

    ! Set up a 3x3-point Gauss rule for 2D
    ! Remark: We will use the 1D rule to create the 2D rule here
    k = 1
    do i = 1, NCUB1D
      do j = 1, NCUB1D
        DcubPts2D(1,k) = DcubPts1D(i)
        DcubPts2D(2,k) = DcubPts1D(j)
        DcubOmega2D(k) = DcubOmega1D(i)*DcubOmega1D(j)
        k = k+1
      end do
    end do


    ! Loop over all elements
    !$omp parallel do default(shared)&
    !$omp private(Da,Dc,DedgeLen,DedgePoints,DedgeWeights,DjacPrep,DjacTrafo,&
    !$omp         DquadPoints,DquadWeights,Dr,Ds,Dvert,bsuccess,ddets,derx,dery,&
    !$omp         dquadArea,dt,dx,dy,i,ipt,j)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Calculate vertice and edge midpoint coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      Dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      Dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      Ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - Dr(2)
      Ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - Dr(1))
      Ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - Dr(2))
      Ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - Dr(1)
      ddets = 1.0_DP / (Ds(1,1)*Ds(2,2) - Ds(1,2)*Ds(2,1))
      Ds(1,1) = ddets*Ds(1,1)
      Ds(1,2) = ddets*Ds(1,2)
      Ds(2,1) = ddets*Ds(2,1)
      Ds(2,2) = ddets*Ds(2,2)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Map 1D cubature points onto the real edges
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Map the 1D cubature points onto the real edges and calculate the
      ! integration weighting factors in this step.
      do j = 1, 4
        ! jacobi determinant of the mapping
        dt = 0.5_DP * sqrt((Dvert(1,mod(j,4)+1)-Dvert(1,j))**2 &
                          +(Dvert(2,mod(j,4)+1)-Dvert(2,j))**2)
        do i = 1, NCUB1D
          DedgePoints(1,i,j) = Dvert(1,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(1,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgePoints(2,i,j) = Dvert(2,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(2,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgeWeights(i,j) = dt * DcubOmega1D(i)
        end do
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Map 2D cubature points onto the real element
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Map the 2D cubature points onto the real element and calculate the
      ! integration weighting factors in this step.
      ! TODO: Replace by Q2-mapping later.
      call trafo_prepJac_quad2D(Dvert, DjacPrep)
      do i = 1, NCUB2D
        call trafo_calcTrafo_quad2d(DjacPrep, DjacTrafo, dt, &
            DcubPts2D(1,i), DcubPts2D(2,i), DquadPoints(1,i), DquadPoints(2,i))
        DquadWeights(i) = dt * DcubOmega2D(i)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 5: Calculate edge lengths and quad area
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the inverse of the edge lengths - we will need them for
      ! scaling later...
      do j = 1, 4
        dt = 0.0_DP
        do i = 1, NCUB1D
          dt = dt + DedgeWeights(i,j)
        end do
        DedgeLen(j) = 1.0_DP / dt
      end do

      ! ...and also calculate the inverse of the element`s area.
      dt = 0.0_DP
      do i = 1, NCUB2D
        dt = dt + DquadWeights(i)
      end do
      dquadArea = 1.0_DP / dt

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 6: Build coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Clear coefficient matrix
      Da = 0.0_DP

      ! Loop over all edges of the quad
      do j = 1, 4

        ! Loop over all cubature points on the current edge
        do i = 1, NCUB1D

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(DedgePoints(1,i,j)-Dr(1)) &
             + Ds(1,2)*(DedgePoints(2,i,j)-Dr(2))
          dy = Ds(2,1)*(DedgePoints(1,i,j)-Dr(1)) &
             + Ds(2,2)*(DedgePoints(2,i,j)-Dr(2))

          ! Integral-Mean over the edges
          ! ----------------------------
          dt = DedgeWeights(i,j) * DedgeLen(j)

          Da(1,j) = Da(1,j) + dt
          Da(2,j) = Da(2,j) + dx*dt
          Da(3,j) = Da(3,j) + dy*dt
          Da(4,j) = Da(4,j) + dx*dy*dt
          Da(5,j) = Da(5,j) + (dx**2 - dy**2)*dt

        end do ! i

      end do ! j

      ! Loop over all 2D cubature points
      do i = 1, NCUB2D

        ! Apply inverse affine trafo to get (x,y)
        dx = Ds(1,1)*(DquadPoints(1,i)-Dr(1)) &
           + Ds(1,2)*(DquadPoints(2,i)-Dr(2))
        dy = Ds(2,1)*(DquadPoints(1,i)-Dr(1)) &
           + Ds(2,2)*(DquadPoints(2,i)-Dr(2))

        ! Integral-Mean over the element
        ! ------------------------------
        dt = DquadWeights(i) * dquadArea *dx * dy

        Da(1,5) = Da(1,5) + dt
        Da(2,5) = Da(2,5) + dx*dt
        Da(3,5) = Da(3,5) + dy*dt
        Da(4,5) = Da(4,5) + dx*dy*dt
        Da(5,5) = Da(5,5) + (dx**2 - dy**2)*dt

      end do ! i

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 7: Invert coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Call the 'direct' inversion routine for 5x5 systems
      call mprim_invert5x5MatrixDirect(Da, Dc, bsuccess)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 8: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))
          dy = Ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))

          ! Evaluate basis functions
          do i = 1, NBAS

            Dbas(i,DER_FUNC2D,ipt,iel) = Dc(i,1) + dx*(Dc(i,2) + dx*Dc(i,5)) &
                                 + dx*dy*Dc(i,4) + dy*(Dc(i,3) - dy*Dc(i,5))

          end do ! i

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 9: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))
          dy = Ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))

          ! Evaluate derivatives
          do i = 1, NBAS

            ! Calculate 'reference' derivatives
            derx = Dc(i,2) + dy*Dc(i,4) + 2.0_DP*dx*Dc(i,5)
            dery = Dc(i,3) + dx*Dc(i,4) - 2.0_DP*dy*Dc(i,5)

            ! Calculate 'real' derivatives
            Dbas(i,DER_DERIV2D_X,ipt,iel) = Ds(1,1)*derx + Ds(2,1)*dery
            Dbas(i,DER_DERIV2D_Y,ipt,iel) = Ds(1,2)*derx + Ds(2,2)*dery

          end do ! i

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_EN30_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The EN30_2D element is specified by four polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x^2 - y^2 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,4:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))) d(t) = kronecker(i,j) * |ej|
  !     <==>
  !     Int_ej (Pi(t)) d(t) = kronecker(i,j) * |ej|
  !   }
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t


  ! Parameter: Number of local basis functions
  integer, parameter :: NBAS = 4

  ! Parameter: Number of cubature points for 1D edge integration
  integer, parameter :: NCUB1D = 2

  ! 1D edge cubature rule point coordinates and weights
  real(DP), dimension(NCUB1D) :: DcubPts1D
  real(DP), dimension(NCUB1D) :: DcubOmega1D

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert

  ! Local mapped 1D cubature point coordinates and integration weights
  real(DP), dimension(NDIM2D, NCUB1D, 4) :: DedgePoints
  real(DP), dimension(NCUB1D, 4) :: DedgeWeights
  real(DP), dimension(4) :: DedgeLen

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr
  real(DP) :: ddets

  ! other local variables
  logical :: bsuccess
  integer :: i,j,iel,ipt
  real(DP), dimension(NBAS,NBAS) :: Da, Dc
  real(DP) :: dx,dy,dt,derx,dery


    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Step 0: Set up 1D cubature rules
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Set up a 2-point Gauss rule for 1D
    DcubPts1D(1) = -sqrt(1.0_DP / 3.0_DP)
    DcubPts1D(2) =  sqrt(1.0_DP / 3.0_DP)
    DcubOmega1D(1) = 1.0_DP
    DcubOmega1D(2) = 1.0_DP

    ! Loop over all elements
    !$omp parallel do default(shared)&
    !$omp private(Da,Dc,DedgeLen,DedgePoints,DedgeWeights,Dr,Ds,Dvert,bsuccess,&
    !$omp         ddets,derx,dery,dt,dx,dy,i,ipt,j)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Fetch vertice coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - dr(2)
      ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - dr(1))
      ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - dr(2))
      ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - dr(1)
      ddets = 1.0_DP / (ds(1,1)*ds(2,2) - ds(1,2)*ds(2,1))
      Ds = ddets * Ds

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Map 1D cubature points onto the real edges
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Map the 1D cubature points onto the real edges and calculate the
      ! integration weighting factors in this step.
      do j = 1, 4
        ! jacobi determinant of the mapping
        dt = 0.5_DP * sqrt((Dvert(1,mod(j,4)+1)-Dvert(1,j))**2 &
                          +(Dvert(2,mod(j,4)+1)-Dvert(2,j))**2)
        do i = 1, NCUB1D
          DedgePoints(1,i,j) = Dvert(1,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(1,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgePoints(2,i,j) = Dvert(2,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(2,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgeWeights(i,j) = dt * DcubOmega1D(i)
        end do
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Calculate edge lengths
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the inverse of the edge lengths - we will need them for
      ! scaling later...
      do j = 1, 4
        dt = 0.0_DP
        do i = 1, NCUB1D
          dt = dt + DedgeWeights(i,j)
        end do
        DedgeLen(j) = 1.0_DP / dt
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 5: Build coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Clear coefficient matrix
      Da = 0.0_DP

      ! Loop over all edges of the quad
      do j = 1, 4

        ! Loop over all cubature points on the current edge
        do i = 1, NCUB1D

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(DedgePoints(1,i,j)-dr(1)) &
             + ds(1,2)*(DedgePoints(2,i,j)-dr(2))
          dy = ds(2,1)*(DedgePoints(1,i,j)-dr(1)) &
             + ds(2,2)*(DedgePoints(2,i,j)-dr(2))

          ! Integral-Mean over the edges
          ! ----------------------------
          dt = DedgeWeights(i,j) * DedgeLen(j)

          ! Evaluate m1(x,y) = 1
          Da(1,j) = Da(1,j) + dt
          ! Evaluate m2(x,y) = x
          Da(2,j) = Da(2,j) + dx*dt
          ! Evaluate m3(x,y) = y
          Da(3,j) = Da(3,j) + dy*dt
          ! Evaluate m4(x,y) = x^2 - y^2
          Da(4,j) = Da(4,j) + (dx**2 - dy**2)*dt

        end do ! i

      end do ! j

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 6: Invert coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Call the 'direct' inversion routine for 4x4 systems
      call mprim_invert4x4MatrixDirect(Da, Dc, bsuccess)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 7: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate basis functions
          do i = 1, NBAS

            Dbas(i,DER_FUNC2D,ipt,iel) = Dc(i,1) + dx*(Dc(i,2) + dx*Dc(i,4)) &
                                                 + dy*(Dc(i,3) - dy*Dc(i,4))
          end do ! i

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 9: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate derivatives
          do i = 1, NBAS

            ! Calculate 'reference' derivatives
            derx = Dc(i,2) + 2.0_DP*dx*Dc(i,4)
            dery = Dc(i,3) - 2.0_DP*dy*Dc(i,4)

            ! Calculate 'real' derivatives
            Dbas(i,DER_DERIV2D_X,ipt,iel) = ds(1,1)*derx + ds(2,1)*dery
            Dbas(i,DER_DERIV2D_Y,ipt,iel) = ds(1,2)*derx + ds(2,2)*dery

          end do ! i

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_EN31_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The EM31_2D element is specified by four polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x^2 - y^2 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,4:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Pi(ej) = kronecker(i,j)
  !   }
  ! }
  !
  ! With:
  ! ej being the midpoint of the j-th local edge of the quadrilateral
  !
  ! Although this element implementation is non-parametric, we can use
  ! exactly the same basis functions as the parametric variant, since
  ! the edge midpoints of the real element are mapped onto the edge
  ! midpoints of the reference element.
  !
  !  P1(x,y) = -(x^2-y^2)/4 - y/2 + 1/4
  !  P2(x,y) =  (x^2-y^2)/4 + x/2 + 1/4
  !  P3(x,y) = -(x^2-y^2)/4 + y/2 + 1/4
  !  P4(x,y) =  (x^2-y^2)/4 - x/2 + 1/4


  ! Parameter: Number of local basis functions
  integer, parameter :: NBAS = 4

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr

  ! other local variables
  integer :: i,iel,ipt
  real(DP), dimension(NDIM2D,NBAS) :: Dgrad
  real(DP) :: dx,dy,dt

    ! Loop over all elements
    !$omp parallel do default(shared)&
    !$omp private(Dgrad,Dr,Ds,Dvert,dt,dx,dy,i,ipt)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Fetch vertice coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - dr(2)
      ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - dr(1))
      ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - dr(2))
      ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - dr(1)
      dt = 1.0_DP / (ds(1,1)*ds(2,2) - ds(1,2)*ds(2,1))
      Ds = dt * Ds

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate basis functions
          dt = (dx**2 - dy**2)
          Dbas(1,DER_FUNC2D,ipt,iel) = 0.25_DP*(-dt - 2.0_DP*dy + 1.0_DP)
          Dbas(2,DER_FUNC2D,ipt,iel) = 0.25_DP*( dt + 2.0_DP*dx + 1.0_DP)
          Dbas(3,DER_FUNC2D,ipt,iel) = 0.25_DP*(-dt + 2.0_DP*dy + 1.0_DP)
          Dbas(4,DER_FUNC2D,ipt,iel) = 0.25_DP*( dt - 2.0_DP*dx + 1.0_DP)

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Calculate 'reference' derivatives
          Dgrad(1,1) = -0.5_DP*dx
          Dgrad(1,2) =  0.5_DP*(dx + 1.0_DP)
          Dgrad(1,3) = -0.5_DP*dx
          Dgrad(1,4) =  0.5_DP*(dx - 1.0_DP)
          Dgrad(2,1) =  0.5_DP*(dy - 1.0_DP)
          Dgrad(2,2) = -0.5_DP*dy
          Dgrad(2,3) =  0.5_DP*(dy + 1.0_DP)
          Dgrad(2,4) = -0.5_DP*dy

          ! Calculate 'real' derivatives
          do i = 1, NBAS
            Dbas(i,DER_DERIV2D_X,ipt,iel) = ds(1,1)*Dgrad(1,i) + ds(2,1)*Dgrad(2,i)
            Dbas(i,DER_DERIV2D_Y,ipt,iel) = ds(1,2)*Dgrad(1,i) + ds(2,2)*Dgrad(2,i)
          end do ! i

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_E032_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The E032_2D element is specified by five polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x^2, y^2 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,5:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))) d(t) = kronecker(i,j) * |ej|
  !     <==>
  !     Int_ej (Pi(t)) d(t) = kronecker(i,j) * |ej|
  !   }
  !   Int_T (Pi(x,y)) d(x,y) = kronecker(i,5) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  ! T being the quadrilateral
  ! |T| being the area of the quadrilateral
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  !  P1(x,y) = ((3*y - 2)*y - 1)*1/4
  !  P2(x,y) = ((3*x + 2)*x - 1)*1/4
  !  P3(x,y) = ((3*y + 2)*y - 1)*1/4
  !  P4(x,y) = ((3*x - 2)*x - 1)*1/4
  !  P5(x,y) = 2 - 3/2*(x^2 + y^2)

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 5

  ! Local variables
  real(DP) :: ddet,dx,dy
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = ((3.0_DP*dy - 2.0_DP)*dy - 1.0_DP)*0.25_DP
          Dbas(2,DER_FUNC2D,i,j) = ((3.0_DP*dx + 2.0_DP)*dx - 1.0_DP)*0.25_DP
          Dbas(3,DER_FUNC2D,i,j) = ((3.0_DP*dy + 2.0_DP)*dy - 1.0_DP)*0.25_DP
          Dbas(4,DER_FUNC2D,i,j) = ((3.0_DP*dx - 2.0_DP)*dx - 1.0_DP)*0.25_DP
          Dbas(5,DER_FUNC2D,i,j) = 2.0_DP - 1.5_DP*(dx*dx + dy*dy)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(DrefDer,i,ddet,dx,dy)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) = 0.0_DP
          DrefDer(2,1) = 1.5_DP*dx + 0.5_DP
          DrefDer(3,1) = 0.0_DP
          DrefDer(4,1) = 1.5_DP*dx - 0.5_DP
          DrefDer(5,1) = -3.0_DP*dx
          ! Y-derivatives
          DrefDer(1,2) = 1.5_DP*dy - 0.5_DP
          DrefDer(2,2) = 0.0_DP
          DrefDer(3,2) = 1.5_DP*dy + 0.5_DP
          DrefDer(4,2) = 0.0_DP
          DrefDer(5,2) = -3.0_DP*dy

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_E050_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The E050_2D element is specified by nine polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2, y^2, x^2*y, x*y^2, x^3*y - x*y^3 }
  !
  ! see:
  ! J.-P. Hennart, J. Jaffre, and J. E. Roberts;
  ! "A Constructive Method for Deriving Finite Elements of Nodal Type";
  ! Numer. Math., 53 (1988), pp. 701-738.
  ! (The basis monomial set above is presented in example 5, pp. 716-720;
  !  see Fig. 8 on page 717)
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,9:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))      ) d(t) = kronecker(i,j  ) * |ej|
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))*L1(t)) d(t) = kronecker(i,j+4) * |ej|
  !   }
  !   Int_T (Pi(x,y)) d(x,y) = kronecker(i,9) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  ! T being the quadrilateral
  ! |T| being the area of the quadrilateral
  ! L1 being the first Legendre-Polynomial:
  ! L1(x) := x
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  !  P1 (x,y) =  3/4*y*( x^2 + y - 1) - 1/4
  !  P2 (x,y) =  3/4*x*(-y^2 + x + 1) - 1/4
  !  P3 (x,y) =  3/4*y*(-x^2 + y + 1) - 1/4
  !  P4 (x,y) =  3/4*x*( y^2 + x - 1) - 1/4
  !  P5 (x,y) = -3/8*x*(y*(y*( 5*y - 6) - 5*x^2 + 2) + 2)
  !  P6 (x,y) = -3/8*y*(x*(x*(-5*x - 6) + 5*y^2 - 2) + 2)
  !  P7 (x,y) = -3/8*x*(y*(y*( 5*y + 6) - 5*x^2 + 2) - 2)
  !  P8 (x,y) = -3/8*y*(x*(x*(-5*x + 6) + 5*y^2 - 2) - 2)
  !  P9 (x,y) = -3/2*(x^2 + y^2) + 2

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 9

  ! Local variables
  real(DP) :: ddet,dx,dy,dx2,dy2
  integer :: i,j
  integer(I32) :: itwist
  real(DP), dimension(4) :: Dtw

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: Q8 = 6.25_DP  ! = 25/4
  real(DP), parameter :: Q9 = 37.5_DP  ! = 75/2
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Get the twist indices for this element.
        itwist = reval%p_ItwistIndex(j)
        Dtw(1) = real(1-iand(int(ishft(itwist, 1)),2),DP)
        Dtw(2) = real(1-iand(int(ishft(itwist, 0)),2),DP)
        Dtw(3) = real(1-iand(int(ishft(itwist,-1)),2),DP)
        Dtw(4) = real(1-iand(int(ishft(itwist,-2)),2),DP)

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Pre-calculate squares
          dx2 = dx*dx
          dy2 = dy*dy

          ! Evaluate basis functions
          Dbas( 1,DER_FUNC2D,i,j) =  Q3*dy*( dx2 + dy - P1) - Q1
          Dbas( 2,DER_FUNC2D,i,j) =  Q3*dx*(-dy2 + dx + P1) - Q1
          Dbas( 3,DER_FUNC2D,i,j) =  Q3*dy*(-dx2 + dy + P1) - Q1
          Dbas( 4,DER_FUNC2D,i,j) =  Q3*dx*( dy2 + dx - P1) - Q1
          Dbas( 5,DER_FUNC2D,i,j) = -Q2*dx*(dy*(dy*( P4*dy - P5) &
                                          - P4*dx2 + P2) + P2)*Dtw(1)
          Dbas( 6,DER_FUNC2D,i,j) = -Q2*dy*(dx*(dx*(-P4*dx - P5) &
                                          + P4*dy2 - P2) + P2)*Dtw(2)
          Dbas( 7,DER_FUNC2D,i,j) = -Q2*dx*(dy*(dy*( P4*dy + P5) &
                                          - P4*dx2 + P2) - P2)*Dtw(3)
          Dbas( 8,DER_FUNC2D,i,j) = -Q2*dy*(dx*(dx*(-P4*dx + P5) &
                                          + P4*dy2 - P2) - P2)*Dtw(4)
          Dbas( 9,DER_FUNC2D,i,j) = -Q4*(dx2 + dy2) + P2

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2,ddet,DrefDer) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Get the twist indices for this element.
        itwist = reval%p_ItwistIndex(j)
        Dtw(1) = real(1-iand(int(ishft(itwist, 1)),2),DP)
        Dtw(2) = real(1-iand(int(ishft(itwist, 0)),2),DP)
        Dtw(3) = real(1-iand(int(ishft(itwist,-1)),2),DP)
        Dtw(4) = real(1-iand(int(ishft(itwist,-2)),2),DP)

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Pre-calculate squares
          dx2 = dx*dx
          dy2 = dy*dy

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer( 1,1) =  Q4*dx*dy
          DrefDer( 2,1) = -Q3*dy2 + Q4*dx + Q3
          DrefDer( 3,1) = -Q4*dx*dy
          DrefDer( 4,1) =  Q3*dy2 + Q4*dx - Q3
          DrefDer( 5,1) = ( dy*(dy*(-Q6*dy + Q5) + Q7*dx2 - Q3) - Q3)*Dtw(1)
          DrefDer( 6,1) = (-Q2*dy*(dx*(-P7*dx - P6) + P4*dy2 - P2))*Dtw(2)
          DrefDer( 7,1) = ( dy*(dy*(-Q6*dy - Q5) + Q7*dx2 - Q3) + Q3)*Dtw(3)
          DrefDer( 8,1) = (-Q2*dy*(dx*(-P7*dx + P6) + P4*dy2 - P2))*Dtw(4)
          DrefDer( 9,1) = -P3*dx
          ! Y-derivatives
          DrefDer( 1,2) =  Q3*dx2 + Q4*dy - Q3
          DrefDer( 2,2) = -Q4*dx*dy
          DrefDer( 3,2) = -Q3*dx2 + Q4*dy + Q3
          DrefDer( 4,2) =  Q4*dx*dy
          DrefDer( 5,2) = (-Q2*dx*(dy*( P7*dy - P6) - P4*dx2 + P2))*Dtw(1)
          DrefDer( 6,2) = ( dx*(dx*( Q6*dx + Q5) - Q7*dy2 + Q3) - Q3)*Dtw(2)
          DrefDer( 7,2) = ( Q2*dx*(dy*(-P7*dy - P6) + P4*dx2 - P2))*Dtw(3)
          DrefDer( 8,2) = ( dx*(dx*( Q6*dx - Q5) - Q7*dy2 + Q3) + Q3)*Dtw(4)
          DrefDer( 9,2) = -P3*dy

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_EB50_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The EB50_2D element is specified by ten polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2, y^2, x^2*y, x*y^2, x^2*y^2, x^3*y - x*y^3 }
  !
  ! see:
  ! J.-P. Hennart, J. Jaffre, and J. E. Roberts;
  ! "A Constructive Method for Deriving Finite Elements of Nodal Type";
  ! Numer. Math., 53 (1988), pp. 701-738.
  ! (The basis monomial set above is an extension of the monomial set
  !  presented in example 5, pp. 716-720; see Fig. 8 on page 717)
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,10:
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))      ) d(t) = kronecker(i,j  ) * |ej|
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))*L1(t)) d(t) = kronecker(i,j+4) * |ej|
  !   }
  !   Int_T (Pi(x,y)            ) d(x,y) = kronecker(i, 9) * |T|
  !   Int_T (Pi(x,y)*L2(x)*L2(y)) d(x,y) = kronecker(i,10) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  ! T being the quadrilateral
  ! |T| being the area of the quadrilateral
  ! L1 and L2 being the first two Legendre-Polynomials:
  ! L1(x) := x
  ! L2(x) := 1/2*(3*x^2 - 1)
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  !  P1 (x,y) =  3/4*y*( x^2 + y - 1) - 1/4
  !  P2 (x,y) =  3/4*x*(-y^2 + x + 1) - 1/4
  !  P3 (x,y) =  3/4*y*(-x^2 + y + 1) - 1/4
  !  P4 (x,y) =  3/4*x*( y^2 + x - 1) - 1/4
  !  P5 (x,y) = -3/8*x*(y*(y*( 5*y - 6) - 5*x^2 + 2) + 2)
  !  P6 (x,y) = -3/8*y*(x*(x*(-5*x - 6) + 5*y^2 - 2) + 2)
  !  P7 (x,y) = -3/8*x*(y*(y*( 5*y + 6) - 5*x^2 + 2) - 2)
  !  P8 (x,y) = -3/8*y*(x*(x*(-5*x + 6) + 5*y^2 - 2) - 2)
  !  P9 (x,y) = -3/2*(x^2 + y^2) + 2
  !  P10(x,y) = 25/4*(3*x^2 - 1)*(3*y^2 - 1)

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 10

  ! Local variables
  real(DP) :: ddet,dx,dy,dx2,dy2
  integer :: i,j
  integer(I32) :: itwist
  real(DP), dimension(4) :: Dtw

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM2D) :: DrefDer

  ! A hand full of parameters to make the code less readable ^_^
  real(DP), parameter :: Q1 = 0.25_DP  ! =  1/4
  real(DP), parameter :: Q2 = 0.375_DP ! =  3/8
  real(DP), parameter :: Q3 = 0.75_DP  ! =  3/4
  real(DP), parameter :: Q4 = 1.5_DP   ! =  3/2
  real(DP), parameter :: Q5 = 2.25_DP  ! =  9/4
  real(DP), parameter :: Q6 = 1.875_DP ! = 15/8
  real(DP), parameter :: Q7 = 5.625_DP ! = 45/8
  real(DP), parameter :: Q8 = 6.25_DP  ! = 25/4
  real(DP), parameter :: Q9 = 37.5_DP  ! = 75/2
  real(DP), parameter :: P1 = 1.0_DP
  real(DP), parameter :: P2 = 2.0_DP
  real(DP), parameter :: P3 = 3.0_DP
  real(DP), parameter :: P4 = 5.0_DP
  real(DP), parameter :: P5 = 6.0_DP
  real(DP), parameter :: P6 = 12.0_DP
  real(DP), parameter :: P7 = 15.0_DP

    ! Calculate function values?
    if(Bder(DER_FUNC2D)) then

      ! Loop through all elements
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Get the twist indices for this element.
        itwist = reval%p_ItwistIndex(j)
        Dtw(1) = real(1-iand(int(ishft(itwist, 1)),2),DP)
        Dtw(2) = real(1-iand(int(ishft(itwist, 0)),2),DP)
        Dtw(3) = real(1-iand(int(ishft(itwist,-1)),2),DP)
        Dtw(4) = real(1-iand(int(ishft(itwist,-2)),2),DP)

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Pre-calculate squares
          dx2 = dx*dx
          dy2 = dy*dy

          ! Evaluate basis functions
          Dbas( 1,DER_FUNC2D,i,j) =  Q3*dy*( dx2 + dy - P1) - Q1
          Dbas( 2,DER_FUNC2D,i,j) =  Q3*dx*(-dy2 + dx + P1) - Q1
          Dbas( 3,DER_FUNC2D,i,j) =  Q3*dy*(-dx2 + dy + P1) - Q1
          Dbas( 4,DER_FUNC2D,i,j) =  Q3*dx*( dy2 + dx - P1) - Q1
          Dbas( 5,DER_FUNC2D,i,j) = -Q2*dx*(dy*(dy*( P4*dy - P5) &
                                          - P4*dx2 + P2) + P2)*Dtw(1)
          Dbas( 6,DER_FUNC2D,i,j) = -Q2*dy*(dx*(dx*(-P4*dx - P5) &
                                          + P4*dy2 - P2) + P2)*Dtw(2)
          Dbas( 7,DER_FUNC2D,i,j) = -Q2*dx*(dy*(dy*( P4*dy + P5) &
                                          - P4*dx2 + P2) - P2)*Dtw(3)
          Dbas( 8,DER_FUNC2D,i,j) = -Q2*dy*(dx*(dx*(-P4*dx + P5) &
                                          + P4*dy2 - P2) - P2)*Dtw(4)
          Dbas( 9,DER_FUNC2D,i,j) = -Q4*(dx2 + dy2) + P2
          Dbas(10,DER_FUNC2D,i,j) =  Q8*(P3*dx2 - P1)*(P3*dy2 - P1)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2,ddet,DrefDer) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Get the twist indices for this element.
        itwist = reval%p_ItwistIndex(j)
        Dtw(1) = real(1-iand(int(ishft(itwist, 1)),2),DP)
        Dtw(2) = real(1-iand(int(ishft(itwist, 0)),2),DP)
        Dtw(3) = real(1-iand(int(ishft(itwist,-1)),2),DP)
        Dtw(4) = real(1-iand(int(ishft(itwist,-2)),2),DP)

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)

          ! Pre-calculate squares
          dx2 = dx*dx
          dy2 = dy*dy

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer( 1,1) =  Q4*dx*dy
          DrefDer( 2,1) = -Q3*dy2 + Q4*dx + Q3
          DrefDer( 3,1) = -Q4*dx*dy
          DrefDer( 4,1) =  Q3*dy2 + Q4*dx - Q3
          DrefDer( 5,1) = ( dy*(dy*(-Q6*dy + Q5) + Q7*dx2 - Q3) - Q3)*Dtw(1)
          DrefDer( 6,1) = (-Q2*dy*(dx*(-P7*dx - P6) + P4*dy2 - P2))*Dtw(2)
          DrefDer( 7,1) = ( dy*(dy*(-Q6*dy - Q5) + Q7*dx2 - Q3) + Q3)*Dtw(3)
          DrefDer( 8,1) = (-Q2*dy*(dx*(-P7*dx + P6) + P4*dy2 - P2))*Dtw(4)
          DrefDer( 9,1) = -P3*dx
          DrefDer(10,1) =  Q9*dx*(P3*dy2 - P1)
          ! Y-derivatives
          DrefDer( 1,2) =  Q3*dx2 + Q4*dy - Q3
          DrefDer( 2,2) = -Q4*dx*dy
          DrefDer( 3,2) = -Q3*dx2 + Q4*dy + Q3
          DrefDer( 4,2) =  Q4*dx*dy
          DrefDer( 5,2) = (-Q2*dx*(dy*( P7*dy - P6) - P4*dx2 + P2))*Dtw(1)
          DrefDer( 6,2) = ( dx*(dx*( Q6*dx + Q5) - Q7*dy2 + Q3) - Q3)*Dtw(2)
          DrefDer( 7,2) = ( Q2*dx*(dy*(-P7*dy - P6) + P4*dx2 - P2))*Dtw(3)
          DrefDer( 8,2) = ( dx*(dx*( Q6*dx - Q5) - Q7*dy2 + Q3) + Q3)*Dtw(4)
          DrefDer( 9,2) = -P3*dy
          DrefDer(10,2) =  Q9*dy*(P3*dx2 - P1)

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 2D quad elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = reval%p_Djac(4,i,j)*ddet
          dy = reval%p_Djac(2,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

          ! Y-derivatives on real element
          dx = -reval%p_Djac(3,i,j)*ddet
          dy = -reval%p_Djac(1,i,j)*ddet
          Dbas(1:NBAS,DER_DERIV2D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                                         - dy*DrefDer(1:NBAS,2)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_EN50_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier, must be EL_EN50_2D.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The EN50_2D element is specified by nine polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2, y^2, x^2*y, x*y^2, x^3*y - x*y^3 }
  !
  ! see:
  ! J.-P. Hennart, J. Jaffre, and J. E. Roberts;
  ! "A Constructive Method for Deriving Finite Elements of Nodal Type";
  ! Numer. Math., 53 (1988), pp. 701-738.
  ! (The basis monomial set above is presented in example 5, pp. 716-720;
  !  see Fig. 8 on page 717)
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,9
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))      ) d(t) = kronecker(i,j  ) * |ej|
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))*L1(t)) d(t) = kronecker(i,j+4) * |ej|
  !   }
  !   Int_T (Pi(x,y)) d(x,y) = kronecker(i, 9) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  ! T being the quadrilateral
  ! |T| being the area of the quadrilateral
  ! L1 being the first Legendre-Polynomial:
  ! L1(x) := x


  ! Parameter: Number of local basis functions
  integer, parameter :: NBAS = 9

  ! Parameter: Number of cubature points for 1D edge integration
  integer, parameter :: NCUB1D = 3
  !integer, parameter :: NCUB1D = 5

  ! Parameter: Number of cubature points for 2D quad integration
  integer, parameter :: NCUB2D = NCUB1D**2

  ! 1D edge cubature rule point coordinates and weights
  real(DP), dimension(NCUB1D) :: DcubPts1D
  real(DP), dimension(NCUB1D) :: DcubOmega1D

  ! 2D quad cubature rule point coordinates and weights
  real(DP), dimension(NDIM2D, NCUB2D) :: DcubPts2D
  real(DP), dimension(NCUB2D) :: DcubOmega2D

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert
  !real(DP), dimension(NDIM2D, 4) :: Dedge

  ! Quad midpoint coordinates
  !real(DP), dimension(NDIM2D) :: Dquad

  ! Local mapped 1D cubature point coordinates and integration weights
  real(DP), dimension(NDIM2D, NCUB1D, 4) :: DedgePoints
  real(DP), dimension(NCUB1D, 4) :: DedgeWeights
  real(DP), dimension(4) :: DedgeLen

  ! Local mapped 2D cubature point coordinates and jacobian determinants
  real(DP), dimension(NDIM2D, NCUB2D) :: DquadPoints
  real(DP), dimension(NCUB2D) :: DquadWeights
  real(DP) :: dquadArea

  ! temporary variables for trafo call (will be removed later)
  real(DP), dimension(8) :: DjacPrep
  real(DP), dimension(4) :: DjacTrafo

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr
  real(DP) :: ddets

  ! other local variables
  logical :: bsuccess
  integer(I32) :: itwist
  integer :: i,j,k,iel, ipt
  real(DP), dimension(NBAS,NBAS) :: Da
  real(DP) :: dx,dy,dt,derx,dery
  real(DP), dimension(4) :: Dtwist


    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Step 0: Set up 1D and 2D cubature rules
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Set up a 3-point Gauss rule for 1D
    DcubPts1D(1) = -sqrt(3.0_DP / 5.0_DP)
    DcubPts1D(2) = 0.0_DP
    DcubPts1D(3) = sqrt(3.0_DP / 5.0_DP)
    DcubOmega1D(1) = 5.0_DP / 9.0_DP
    DcubOmega1D(2) = 8.0_DP / 9.0_DP
    DcubOmega1D(3) = 5.0_DP / 9.0_DP

!    ! !!! DEBUG: 5-point Gauss rule !!!
!    dt = 2.0_DP*sqrt(10.0_DP / 7.0_DP)
!    DcubPts1D(1) = -sqrt(5.0_DP + dt) / 3.0_DP
!    DcubPts1D(2) = -sqrt(5.0_DP - dt) / 3.0_DP
!    DcubPts1D(3) = 0.0_DP
!    DcubPts1D(4) =  sqrt(5.0_DP - dt) / 3.0_DP
!    DcubPts1D(5) =  sqrt(5.0_DP + dt) / 3.0_DP
!    dt = 13.0_DP*sqrt(70.0_DP)
!    DcubOmega1D(1) = (322.0_DP - dt) / 900.0_DP
!    DcubOmega1D(2) = (322.0_DP + dt) / 900.0_DP
!    DcubOmega1D(3) = 128.0_DP / 225.0_DP
!    DcubOmega1D(4) = (322.0_DP + dt) / 900.0_DP
!    DcubOmega1D(5) = (322.0_DP - dt) / 900.0_DP

    ! Set up a 3x3-point Gauss rule for 2D
    ! Remark: We will use the 1D rule to create the 2D rule here
    k = 1
    do i = 1, NCUB1D
      do j = 1, NCUB1D
        DcubPts2D(1,k) = DcubPts1D(i)
        DcubPts2D(2,k) = DcubPts1D(j)
        DcubOmega2D(k) = DcubOmega1D(i)*DcubOmega1D(j)
        k = k+1
      end do
    end do

    ! Loop over all elements
    !$omp parallel do default(shared)&
    !$omp private(Da,DedgeLen,DedgePoints,DedgeWeights,DjacPrep,DjacTrafo,&
    !$omp         DquadPoints,DquadWeights,Dr,Ds,Dtwist,Dvert,bsuccess,ddets,&
    !$omp         derx,dery,dquadArea,dt,dx,dy,i,ipt,itwist,j)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Calculate vertice and edge midpoint coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! Calculate edge-midpoint coordinates.
      ! Remark: These may be replaced by 'other' edge midpoints to achieve an
      ! 'iso-parametric' behaviour in a later implementation. (todo)
      !Dedge(1:2,1) = 0.5_DP * (Dvert(1:2,1) + Dvert(1:2,2))
      !Dedge(1:2,2) = 0.5_DP * (Dvert(1:2,2) + Dvert(1:2,3))
      !Dedge(1:2,3) = 0.5_DP * (Dvert(1:2,3) + Dvert(1:2,4))
      !Dedge(1:2,4) = 0.5_DP * (Dvert(1:2,4) + Dvert(1:2,1))

      ! Calculate quad-midpoint coordinates.
      ! Remark: This may be replaced by 'another' quad midpoint to achieve an
      ! 'iso-parametric' behaviour in a later implementation. (todo)
      !Dquad(1:2) = 0.25_DP * (Dvert(1:2,1) + Dvert(1:2,2) &
      !                      + Dvert(1:2,3) + Dvert(1:2,4))

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      Dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      Dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      Ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - Dr(2)
      Ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - Dr(1))
      Ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - Dr(2))
      Ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - Dr(1)
      ddets = 1.0_DP / (Ds(1,1)*Ds(2,2) - Ds(1,2)*Ds(2,1))
      Ds(1,1) = ddets*Ds(1,1)
      Ds(1,2) = ddets*Ds(1,2)
      Ds(2,1) = ddets*Ds(2,1)
      Ds(2,2) = ddets*Ds(2,2)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Map 1D cubature points onto the real edges
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Map the 1D cubature points onto the real edges and calculate the
      ! integration weighting factors in this step.
      ! TODO: Replace by P2-mapping later.
      do j = 1, 4
        ! jacobi determinant of the mapping
        dt = 0.5_DP * sqrt((Dvert(1,mod(j,4)+1)-Dvert(1,j))**2 &
                          +(Dvert(2,mod(j,4)+1)-Dvert(2,j))**2)
        do i = 1, NCUB1D
          DedgePoints(1,i,j) = Dvert(1,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(1,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgePoints(2,i,j) = Dvert(2,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(2,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgeWeights(i,j) = dt * DcubOmega1D(i)
        end do
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Map 2D cubature points onto the real element
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Map the 2D cubature points onto the real element and calculate the
      ! integration weighting factors in this step.
      ! TODO: Replace by Q2-mapping later.
      call trafo_prepJac_quad2D(Dvert, DjacPrep)
      do i = 1, NCUB2D
        call trafo_calcTrafo_quad2d(DjacPrep, DjacTrafo, dt, &
            DcubPts2D(1,i), DcubPts2D(2,i), DquadPoints(1,i), DquadPoints(2,i))
        DquadWeights(i) = dt * DcubOmega2D(i)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 5: Calculate edge lengths and quad area
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the inverse of the edge lengths - we will need them for
      ! scaling later...
      do j = 1, 4
        dt = 0.0_DP
        do i = 1, NCUB1D
          dt = dt + DedgeWeights(i,j)
        end do
        DedgeLen(j) = 1.0_DP / dt
      end do

      ! ...and also calculate the inverse of the element`s area.
      dt = 0.0_DP
      do i = 1, NCUB2D
        dt = dt + DquadWeights(i)
      end do
      dquadArea = 1.0_DP / dt

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 6: Build coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Prepare twist indices
      itwist = reval%p_ItwistIndex(iel)
      Dtwist(1) = real(1-iand(int(ishft(itwist, 1)),2),DP)
      Dtwist(2) = real(1-iand(int(ishft(itwist, 0)),2),DP)
      Dtwist(3) = real(1-iand(int(ishft(itwist,-1)),2),DP)
      Dtwist(4) = real(1-iand(int(ishft(itwist,-2)),2),DP)

      ! Clear coefficient matrix
      Da = 0.0_DP

      ! Loop over all edges of the quad
      do j = 1, 4

        ! Loop over all cubature points on the current edge
        do i = 1, NCUB1D

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(DedgePoints(1,i,j)-Dr(1)) &
             + Ds(1,2)*(DedgePoints(2,i,j)-Dr(2))
          dy = Ds(2,1)*(DedgePoints(1,i,j)-Dr(1)) &
             + Ds(2,2)*(DedgePoints(2,i,j)-Dr(2))

          ! Integral-Mean over the edges
          ! ----------------------------
          dt = DedgeWeights(i,j) * DedgeLen(j)

          Da(1,j) = Da(1,j) + dt
          Da(2,j) = Da(2,j) + dx*dt
          Da(3,j) = Da(3,j) + dy*dt
          Da(4,j) = Da(4,j) + dx*dy*dt
          Da(5,j) = Da(5,j) + dx**2*dt
          Da(6,j) = Da(6,j) + dy**2*dt
          Da(7,j) = Da(7,j) + dx**2*dy*dt
          Da(8,j) = Da(8,j) + dx*dy**2*dt
          Da(9,j) = Da(9,j) + (dx**3*dy - dx*dy**3)*dt

          ! Legendre-Weighted Integral-Mean over the edges
          ! ----------------------------------------------
          dt = DedgeWeights(i,j) * Dtwist(j) * DcubPts1D(i) * DedgeLen(j)

          Da(1,j+4) = Da(1,j+4) + dt
          Da(2,j+4) = Da(2,j+4) + dx*dt
          Da(3,j+4) = Da(3,j+4) + dy*dt
          Da(4,j+4) = Da(4,j+4) + dx*dy*dt
          Da(5,j+4) = Da(5,j+4) + dx**2*dt
          Da(6,j+4) = Da(6,j+4) + dy**2*dt
          Da(7,j+4) = Da(7,j+4) + dx**2*dy*dt
          Da(8,j+4) = Da(8,j+4) + dx*dy**2*dt
          Da(9,j+4) = Da(9,j+4) + (dx**3*dy - dx*dy**3)*dt

        end do ! i

      end do ! j

      ! Loop over all 2D cubature points
      do i = 1, NCUB2D

        ! Apply inverse affine trafo to get (x,y)
        dx = Ds(1,1)*(DquadPoints(1,i)-Dr(1)) &
           + Ds(1,2)*(DquadPoints(2,i)-Dr(2))
        dy = Ds(2,1)*(DquadPoints(1,i)-Dr(1)) &
           + Ds(2,2)*(DquadPoints(2,i)-Dr(2))

        ! Integral-Mean over the element
        ! ------------------------------
        dt = DquadWeights(i) * dquadArea

        Da(1,9) = Da(1,9) + dt
        Da(2,9) = Da(2,9) + dx*dt
        Da(3,9) = Da(3,9) + dy*dt
        Da(4,9) = Da(4,9) + dx*dy*dt
        Da(5,9) = Da(5,9) + dx**2*dt
        Da(6,9) = Da(6,9) + dy**2*dt
        Da(7,9) = Da(7,9) + dx**2*dy*dt
        Da(8,9) = Da(8,9) + dx*dy**2*dt
        Da(9,9) = Da(9,9) + (dx**3*dy - dx*dy**3)*dt

      end do ! i

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 7: Invert coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call mprim_invertMatrixPivot(Da, NBAS, bsuccess)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 8: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))
          dy = Ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))

          ! Evaluate basis functions
          do i = 1, NBAS

            Dbas(i,DER_FUNC2D,ipt,iel) = Da(i,1) + dx*dy*Da(i,4) &
              + dx*(Da(i,2) + dx*(Da(i,5) + dy*(Da(i,7) + dx*Da(i,9)))) &
              + dy*(Da(i,3) + dy*(Da(i,6) + dx*(Da(i,8) - dy*Da(i,9))))

          end do ! i

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 9: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))
          dy = Ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))

          ! Evaluate derivatives
          do i = 1, NBAS

            ! Calculate 'reference' derivatives
            derx = Da(i,2) + dy*(Da(i,4) + dy*(Da(i,8) - dy*Da(i,9))) &
                 + 2.0_DP*dx*(Da(i,5) + dy*(Da(i,7) + 1.5_DP*dx*Da(i,9)))
            dery = Da(i,3) + dx*(Da(i,4) + dx*(Da(i,7) + dx*Da(i,9))) &
                 + 2.0_DP*dy*(Da(i,6) + dx*(Da(i,8) - 1.5_DP*dy*Da(i,9)))

            ! Calculate 'real' derivatives
            Dbas(i,DER_DERIV2D_X,ipt,iel) = Ds(1,1)*derx + Ds(2,1)*dery
            Dbas(i,DER_DERIV2D_Y,ipt,iel) = Ds(1,2)*derx + Ds(2,2)*dery

          end do ! i

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_EN51_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier, must be EL_EN51_2D.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The EN51_2D element is specified by fifteen polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^2*y^2,
  !   x^3*y^2, x^2*y^3, x^3*y - x*y^3, x^4*y^2 - x^2*y^4 }
  !
  ! see:
  ! J.-P. Hennart, J. Jaffre, and J. E. Roberts;
  ! "A Constructive Method for Deriving Finite Elements of Nodal Type";
  ! Numer. Math., 53 (1988), pp. 701-738.
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,15
  ! {
  !   For all j = 1,...,4:
  !   {
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))      ) d(t) = kronecker(i,j  ) * |ej|
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))*L1(t)) d(t) = kronecker(i,j+4) * |ej|
  !     Int_[-1,1] (|DEj(t)|*Pi(Ej(t))*L2(t)) d(t) = kronecker(i,j+8) * |ej|
  !   }
  !   Int_T (Pi(x,y)      ) d(x,y) = kronecker(i, 13) * |T|
  !   Int_T (Pi(x,y)*L1(x)) d(x,y) = kronecker(i, 14) * |T|
  !   Int_T (Pi(x,y)*L1(y)) d(x,y) = kronecker(i, 15) * |T|
  ! }
  !
  ! With:
  ! ej being the j-th local edge of the quadrilateral
  ! |ej| being the length of the edge ej
  ! Ej: [-1,1] -> ej being the parametrisation of the edge ej
  ! |DEj(t)| being the determinant of the Jacobi-Matrix of Ej in the point t
  ! T being the quadrilateral
  ! |T| being the area of the quadrilateral
  ! L1 and L2 being the first two Legendre-Polynomials:
  ! L1(x) := x
  ! L2(x) := (3*x^2 - 1) / 2


  ! Parameter: Number of local basis functions
  integer, parameter :: NBAS = 15

  ! Parameter: Number of cubature points for 1D edge integration
  !integer, parameter :: NCUB1D = 3
  integer, parameter :: NCUB1D = 5

  ! Parameter: Number of cubature points for 2D quad integration
  integer, parameter :: NCUB2D = NCUB1D**2

  ! 1D edge cubature rule point coordinates and weights
  real(DP), dimension(NCUB1D) :: DcubPts1D
  real(DP), dimension(NCUB1D) :: DcubOmega1D

  ! 2D quad cubature rule point coordinates and weights
  real(DP), dimension(NDIM2D, NCUB2D) :: DcubPts2D
  real(DP), dimension(NCUB2D) :: DcubOmega2D

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert
  !real(DP), dimension(NDIM2D, 4) :: Dedge

  ! Quad midpoint coordinates
  !real(DP), dimension(NDIM2D) :: Dquad

  ! Local mapped 1D cubature point coordinates and integration weights
  real(DP), dimension(NDIM2D, NCUB1D, 4) :: DedgePoints
  real(DP), dimension(NCUB1D, 4) :: DedgeWeights
  real(DP), dimension(4) :: DedgeLen

  ! Local mapped 2D cubature point coordinates and jacobian determinants
  real(DP), dimension(NDIM2D, NCUB2D) :: DquadPoints
  real(DP), dimension(NCUB2D) :: DquadWeights
  real(DP) :: dquadArea

  ! temporary variables for trafo call (will be removed later)
  real(DP), dimension(8) :: DjacPrep
  real(DP), dimension(4) :: DjacTrafo

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr
  real(DP) :: ddets

  ! other local variables
  logical :: bsuccess
  integer(I32) :: itwist
  integer :: i,j,k,iel, ipt
  real(DP), dimension(NBAS,NBAS) :: Da
  real(DP) :: dx,dy,dt,derx,dery
  real(DP), dimension(4) :: Dtwist


    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Step 0: Set up 1D and 2D cubature rules
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!    ! Set up a 3-point Gauss rule for 1D
!    DcubPts1D(1) = -sqrt(3.0_DP / 5.0_DP)
!    DcubPts1D(2) = 0.0_DP
!    DcubPts1D(3) = sqrt(3.0_DP / 5.0_DP)
!    DcubOmega1D(1) = 5.0_DP / 9.0_DP
!    DcubOmega1D(2) = 8.0_DP / 9.0_DP
!    DcubOmega1D(3) = 5.0_DP / 9.0_DP

    ! !!! DEBUG: 5-point Gauss rule !!!
    dt = 2.0_DP*sqrt(10.0_DP / 7.0_DP)
    DcubPts1D(1) = -sqrt(5.0_DP + dt) / 3.0_DP
    DcubPts1D(2) = -sqrt(5.0_DP - dt) / 3.0_DP
    DcubPts1D(3) = 0.0_DP
    DcubPts1D(4) =  sqrt(5.0_DP - dt) / 3.0_DP
    DcubPts1D(5) =  sqrt(5.0_DP + dt) / 3.0_DP
    dt = 13.0_DP*sqrt(70.0_DP)
    DcubOmega1D(1) = (322.0_DP - dt) / 900.0_DP
    DcubOmega1D(2) = (322.0_DP + dt) / 900.0_DP
    DcubOmega1D(3) = 128.0_DP / 225.0_DP
    DcubOmega1D(4) = (322.0_DP + dt) / 900.0_DP
    DcubOmega1D(5) = (322.0_DP - dt) / 900.0_DP

    ! Set up a 3x3-point Gauss rule for 2D
    ! Remark: We will use the 1D rule to create the 2D rule here
    k = 1
    do i = 1, NCUB1D
      do j = 1, NCUB1D
        DcubPts2D(1,k) = DcubPts1D(i)
        DcubPts2D(2,k) = DcubPts1D(j)
        DcubOmega2D(k) = DcubOmega1D(i)*DcubOmega1D(j)
        k = k+1
      end do
    end do

    ! Loop over all elements
    !$omp parallel do default(shared)&
    !$omp private(Da,DedgeLen,DedgePoints,DedgeWeights,DjacPrep,DjacTrafo,&
    !$omp         DquadPoints,DquadWeights,Dr,Ds,Dtwist,Dvert,bsuccess,ddets,&
    !$omp         derx,dery,dquadArea,dt,dx,dy,i,ipt,itwist,j)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Calculate vertice and edge midpoint coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! Calculate edge-midpoint coordinates.
      ! Remark: These may be replaced by 'other' edge midpoints to achieve an
      ! 'iso-parametric' behaviour in a later implementation. (todo)
      !Dedge(1:2,1) = 0.5_DP * (Dvert(1:2,1) + Dvert(1:2,2))
      !Dedge(1:2,2) = 0.5_DP * (Dvert(1:2,2) + Dvert(1:2,3))
      !Dedge(1:2,3) = 0.5_DP * (Dvert(1:2,3) + Dvert(1:2,4))
      !Dedge(1:2,4) = 0.5_DP * (Dvert(1:2,4) + Dvert(1:2,1))

      ! Calculate quad-midpoint coordinates.
      ! Remark: This may be replaced by 'another' quad midpoint to achieve an
      ! 'iso-parametric' behaviour in a later implementation. (todo)
      !Dquad(1:2) = 0.25_DP * (Dvert(1:2,1) + Dvert(1:2,2) &
      !                      + Dvert(1:2,3) + Dvert(1:2,4))

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      Dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      Dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      Ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - Dr(2)
      Ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - Dr(1))
      Ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - Dr(2))
      Ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - Dr(1)
      ddets = 1.0_DP / (Ds(1,1)*Ds(2,2) - Ds(1,2)*Ds(2,1))
      Ds(1,1) = ddets*Ds(1,1)
      Ds(1,2) = ddets*Ds(1,2)
      Ds(2,1) = ddets*Ds(2,1)
      Ds(2,2) = ddets*Ds(2,2)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Map 1D cubature points onto the real edges
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Map the 1D cubature points onto the real edges and calculate the
      ! integration weighting factors in this step.
      do j = 1, 4
        ! jacobi determinant of the mapping
        dt = 0.5_DP * sqrt((Dvert(1,mod(j,4)+1)-Dvert(1,j))**2 &
                          +(Dvert(2,mod(j,4)+1)-Dvert(2,j))**2)
        do i = 1, NCUB1D
          DedgePoints(1,i,j) = Dvert(1,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(1,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgePoints(2,i,j) = Dvert(2,j)*0.5_DP*(1.0_DP - DcubPts1D(i)) &
                    + Dvert(2,mod(j,4)+1)*0.5_DP*(1.0_DP + DcubPts1D(i))
          DedgeWeights(i,j) = dt * DcubOmega1D(i)
        end do
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Map 2D cubature points onto the real element
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Map the 2D cubature points onto the real element and calculate the
      ! integration weighting factors in this step.
      call trafo_prepJac_quad2D(Dvert, DjacPrep)
      do i = 1, NCUB2D
        call trafo_calcTrafo_quad2d(DjacPrep, DjacTrafo, dt, &
            DcubPts2D(1,i), DcubPts2D(2,i), DquadPoints(1,i), DquadPoints(2,i))
        DquadWeights(i) = dt * DcubOmega2D(i)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 5: Calculate edge lengths and quad area
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the inverse of the edge lengths - we will need them for
      ! scaling later...
      do j = 1, 4
        dt = 0.0_DP
        do i = 1, NCUB1D
          dt = dt + DedgeWeights(i,j)
        end do
        DedgeLen(j) = 1.0_DP / dt
      end do

      ! ...and also calculate the inverse of the element`s area.
      dt = 0.0_DP
      do i = 1, NCUB2D
        dt = dt + DquadWeights(i)
      end do
      dquadArea = 1.0_DP / dt

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 6: Build coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Prepare twist indices
      itwist = reval%p_ItwistIndex(iel)
      Dtwist(1) = real(1-iand(int(ishft(itwist, 1)),2),DP)
      Dtwist(2) = real(1-iand(int(ishft(itwist, 0)),2),DP)
      Dtwist(3) = real(1-iand(int(ishft(itwist,-1)),2),DP)
      Dtwist(4) = real(1-iand(int(ishft(itwist,-2)),2),DP)

      ! Clear coefficient matrix
      Da = 0.0_DP

      ! Loop over all edges of the quad
      do j = 1, 4

        ! Loop over all cubature points on the current edge
        do i = 1, NCUB1D

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(DedgePoints(1,i,j)-Dr(1)) &
             + Ds(1,2)*(DedgePoints(2,i,j)-Dr(2))
          dy = Ds(2,1)*(DedgePoints(1,i,j)-Dr(1)) &
             + Ds(2,2)*(DedgePoints(2,i,j)-Dr(2))

          ! Integral-Mean over the edges
          ! ----------------------------
          dt = DedgeWeights(i,j) * DedgeLen(j)

          Da( 1,j) = Da( 1,j) + dt
          Da( 2,j) = Da( 2,j) + dx*dt
          Da( 3,j) = Da( 3,j) + dy*dt
          Da( 4,j) = Da( 4,j) + dx**2*dt
          Da( 5,j) = Da( 5,j) + dx*dy*dt
          Da( 6,j) = Da( 6,j) + dy**2*dt
          Da( 7,j) = Da( 7,j) + dx**3*dt
          Da( 8,j) = Da( 8,j) + dx**2*dy*dt
          Da( 9,j) = Da( 9,j) + dx*dy**2*dt
          Da(10,j) = Da(10,j) + dy**3*dt
          Da(11,j) = Da(11,j) + dx**2*dy**2*dt
          Da(12,j) = Da(12,j) + dx**3*dy**2*dt
          Da(13,j) = Da(13,j) + dx**2*dy**3*dt
          Da(14,j) = Da(14,j) + (dx**3*dy - dx*dy**3)*dt
          Da(15,j) = Da(15,j) + (dx**4*dy**2 - dx**2*dy**4)*dt

          ! Legendre-Weighted Integral-Mean over the edges
          ! ----------------------------------------------
          dt = DedgeWeights(i,j) * DedgeLen(j) * Dtwist(j) * DcubPts1D(i)

          Da( 1,j+4) = Da( 1,j+4) + dt
          Da( 2,j+4) = Da( 2,j+4) + dx*dt
          Da( 3,j+4) = Da( 3,j+4) + dy*dt
          Da( 4,j+4) = Da( 4,j+4) + dx**2*dt
          Da( 5,j+4) = Da( 5,j+4) + dx*dy*dt
          Da( 6,j+4) = Da( 6,j+4) + dy**2*dt
          Da( 7,j+4) = Da( 7,j+4) + dx**3*dt
          Da( 8,j+4) = Da( 8,j+4) + dx**2*dy*dt
          Da( 9,j+4) = Da( 9,j+4) + dx*dy**2*dt
          Da(10,j+4) = Da(10,j+4) + dy**3*dt
          Da(11,j+4) = Da(11,j+4) + dx**2*dy**2*dt
          Da(12,j+4) = Da(12,j+4) + dx**3*dy**2*dt
          Da(13,j+4) = Da(13,j+4) + dx**2*dy**3*dt
          Da(14,j+4) = Da(14,j+4) + (dx**3*dy - dx*dy**3)*dt
          Da(15,j+4) = Da(15,j+4) + (dx**4*dy**2 - dx**2*dy**4)*dt

          ! Legendre-Weighted Integral-Mean over the edges
          ! ----------------------------------------------
          dt = DedgeWeights(i,j) * DedgeLen(j) * (3.0_DP * DcubPts1D(i)**2 - 1.0_DP) * 0.5_DP

          Da( 1,j+8) = Da( 1,j+8) + dt
          Da( 2,j+8) = Da( 2,j+8) + dx*dt
          Da( 3,j+8) = Da( 3,j+8) + dy*dt
          Da( 4,j+8) = Da( 4,j+8) + dx**2*dt
          Da( 5,j+8) = Da( 5,j+8) + dx*dy*dt
          Da( 6,j+8) = Da( 6,j+8) + dy**2*dt
          Da( 7,j+8) = Da( 7,j+8) + dx**3*dt
          Da( 8,j+8) = Da( 8,j+8) + dx**2*dy*dt
          Da( 9,j+8) = Da( 9,j+8) + dx*dy**2*dt
          Da(10,j+8) = Da(10,j+8) + dy**3*dt
          Da(11,j+8) = Da(11,j+8) + dx**2*dy**2*dt
          Da(12,j+8) = Da(12,j+8) + dx**3*dy**2*dt
          Da(13,j+8) = Da(13,j+8) + dx**2*dy**3*dt
          Da(14,j+8) = Da(14,j+8) + (dx**3*dy - dx*dy**3)*dt
          Da(15,j+8) = Da(15,j+8) + (dx**4*dy**2 - dx**2*dy**4)*dt

        end do ! i

      end do ! j

      ! Loop over all 2D cubature points
      do i = 1, NCUB2D

        ! Apply inverse affine trafo to get (x,y)
        dx = Ds(1,1)*(DquadPoints(1,i)-Dr(1)) &
           + Ds(1,2)*(DquadPoints(2,i)-Dr(2))
        dy = Ds(2,1)*(DquadPoints(1,i)-Dr(1)) &
           + Ds(2,2)*(DquadPoints(2,i)-Dr(2))

        ! Integral-Mean over the element
        ! ------------------------------
        dt = DquadWeights(i) * dquadArea

        Da( 1,13) = Da( 1,13) + dt
        Da( 2,13) = Da( 2,13) + dx*dt
        Da( 3,13) = Da( 3,13) + dy*dt
        Da( 4,13) = Da( 4,13) + dx**2*dt
        Da( 5,13) = Da( 5,13) + dx*dy*dt
        Da( 6,13) = Da( 6,13) + dy**2*dt
        Da( 7,13) = Da( 7,13) + dx**3*dt
        Da( 8,13) = Da( 8,13) + dx**2*dy*dt
        Da( 9,13) = Da( 9,13) + dx*dy**2*dt
        Da(10,13) = Da(10,13) + dy**3*dt
        Da(11,13) = Da(11,13) + dx**2*dy**2*dt
        Da(12,13) = Da(12,13) + dx**3*dy**2*dt
        Da(13,13) = Da(13,13) + dx**2*dy**3*dt
        Da(14,13) = Da(14,13) + (dx**3*dy - dx*dy**3)*dt
        Da(15,13) = Da(15,13) + (dx**4*dy**2 - dx**2*dy**4)*dt

        ! Legendre-Weighted Integral-Mean over the element
        ! ------------------------------
        dt = DquadWeights(i) * dquadArea * dx

        Da( 1,14) = Da( 1,14) + dt
        Da( 2,14) = Da( 2,14) + dx*dt
        Da( 3,14) = Da( 3,14) + dy*dt
        Da( 4,14) = Da( 4,14) + dx**2*dt
        Da( 5,14) = Da( 5,14) + dx*dy*dt
        Da( 6,14) = Da( 6,14) + dy**2*dt
        Da( 7,14) = Da( 7,14) + dx**3*dt
        Da( 8,14) = Da( 8,14) + dx**2*dy*dt
        Da( 9,14) = Da( 9,14) + dx*dy**2*dt
        Da(10,14) = Da(10,14) + dy**3*dt
        Da(11,14) = Da(11,14) + dx**2*dy**2*dt
        Da(12,14) = Da(12,14) + dx**3*dy**2*dt
        Da(13,14) = Da(13,14) + dx**2*dy**3*dt
        Da(14,14) = Da(14,14) + (dx**3*dy - dx*dy**3)*dt
        Da(15,14) = Da(15,14) + (dx**4*dy**2 - dx**2*dy**4)*dt

        ! Legendre-Weighted Integral-Mean over the element
        ! ------------------------------
        dt = DquadWeights(i) * dquadArea * dy

        Da( 1,15) = Da( 1,15) + dt
        Da( 2,15) = Da( 2,15) + dx*dt
        Da( 3,15) = Da( 3,15) + dy*dt
        Da( 4,15) = Da( 4,15) + dx**2*dt
        Da( 5,15) = Da( 5,15) + dx*dy*dt
        Da( 6,15) = Da( 6,15) + dy**2*dt
        Da( 7,15) = Da( 7,15) + dx**3*dt
        Da( 8,15) = Da( 8,15) + dx**2*dy*dt
        Da( 9,15) = Da( 9,15) + dx*dy**2*dt
        Da(10,15) = Da(10,15) + dy**3*dt
        Da(11,15) = Da(11,15) + dx**2*dy**2*dt
        Da(12,15) = Da(12,15) + dx**3*dy**2*dt
        Da(13,15) = Da(13,15) + dx**2*dy**3*dt
        Da(14,15) = Da(14,15) + (dx**3*dy - dx*dy**3)*dt
        Da(15,15) = Da(15,15) + (dx**4*dy**2 - dx**2*dy**4)*dt

      end do ! i

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 7: Invert coefficient matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call mprim_invertMatrixPivot(Da, NBAS, bsuccess)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 8: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))
          dy = Ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))

          ! Evaluate basis functions
          do i = 1, NBAS

            Dbas(i,DER_FUNC2D,ipt,iel) = Da(i, 1) &
             + Da(i, 2) * dx &
             + Da(i, 3) * dy &
             + Da(i, 4) * dx**2 &
             + Da(i, 5) * dx*dy &
             + Da(i, 6) * dy**2 &
             + Da(i, 7) * dx**3 &
             + Da(i, 8) * dx**2*dy &
             + Da(i, 9) * dx*dy**2 &
             + Da(i,10) * dy**3 &
             + Da(i,11) * dx**2*dy**2 &
             + Da(i,12) * dx**3*dy**2 &
             + Da(i,13) * dx**2*dy**3 &
             + Da(i,14) * (dx**3*dy - dx*dy**3) &
             + Da(i,15) * (dx**4*dy**2 - dx**2*dy**4)

          end do ! i

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 9: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = Ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))
          dy = Ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-Dr(1)) &
             + Ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-Dr(2))

          ! Evaluate derivatives
          do i = 1, NBAS

            ! Calculate 'reference' derivatives
            derx = Da(i, 2) &
              + Da(i, 4) * dx * 2.0_DP &
              + Da(i, 5) * dy &
              + Da(i, 7) * dx**2 * 3.0_DP &
              + Da(i, 8) * dx*dy * 2.0_DP &
              + Da(i, 9) * dy**2 &
              + Da(i,11) * dx*dy**2 * 2.0_DP &
              + Da(i,12) * dx**2*dy**2 * 3.0_DP &
              + Da(i,13) * dx*dy**3 * 2.0_DP &
              + Da(i,14) * (dx**2*dy*3.0_DP - dy**3) &
              + Da(i,15) * (dx**3*dy**2*4.0_DP - dx*dy**4*2.0_DP)

            dery = Da(i, 3) &
              + Da(i, 5) * dx &
              + Da(i, 6) * dy * 2.0_DP &
              + Da(i, 8) * dx**2 &
              + Da(i, 9) * dx*dy * 2.0_DP &
              + Da(i,10) * dy**2 * 3.0_DP &
              + Da(i,11) * dx**2*dy * 2.0_DP &
              + Da(i,12) * dx**3*dy * 2.0_DP &
              + Da(i,13) * dx**2*dy**2 * 3.0_DP &
              + Da(i,14) * (dx**3 - dx*dy**2 * 3.0_DP) &
              + Da(i,15) * (dx**4*dy * 2.0_DP - dx**2*dy**3 * 4.0_DP)

            ! Calculate 'real' derivatives
            Dbas(i,DER_DERIV2D_X,ipt,iel) = Ds(1,1)*derx + Ds(2,1)*dery
            Dbas(i,DER_DERIV2D_Y,ipt,iel) = Ds(1,2)*derx + Ds(2,2)*dery

          end do ! i

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_DCQP1_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The DCQP1_2D element is specified by three polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y }
  !
  ! As the QP1 element is discontinous, the basis polynomials do not have to
  ! fulfill any special conditions - they are simply defined as:
  !
  !  P1 (x,y,z) = 1
  !  P2 (x,y,z) = x
  !  P3 (x,y,z) = y

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr
  real(DP) :: ddets

  ! other local variables
  integer :: iel,ipt
  real(DP) :: dx,dy

    ! Loop over all elements
    !$omp parallel do default(shared)&
    !$omp private(Dr,Ds,Dvert,ddets,dx,dy,ipt)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Fetch vertice coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - dr(2)
      ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - dr(1))
      ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - dr(2))
      ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - dr(1)
      ddets = 1.0_DP / (ds(1,1)*ds(2,2) - ds(1,2)*ds(2,1))
      Ds = ddets * Ds

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,ipt,iel) = 1.0_DP
          Dbas(2,DER_FUNC2D,ipt,iel) = dx
          DBas(3,DER_FUNC2D,ipt,iel) = dy

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Evaluate derivatives
          Dbas(1,DER_DERIV2D_X,ipt,iel) = 0.0_DP
          Dbas(1,DER_DERIV2D_Y,ipt,iel) = 0.0_DP
          Dbas(2,DER_DERIV2D_X,ipt,iel) = ds(1,1)
          Dbas(2,DER_DERIV2D_Y,ipt,iel) = ds(1,2)
          Dbas(3,DER_DERIV2D_X,ipt,iel) = ds(2,1)
          Dbas(3,DER_DERIV2D_Y,ipt,iel) = ds(2,2)

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_DCQP2_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The DCQP2_2D element is specified by six polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x^2, y^2, x*y }
  !
  ! As the DCQP2_2D element is discontinous, the basis polynomials do not have to
  ! fulfill any special conditions - they are simply defined as:
  !
  !  P1 (x,y,z) = 1
  !  P2 (x,y,z) = x
  !  P3 (x,y,z) = y
  !  P4 (x,y,z) = x^2
  !  P5 (x,y,z) = y^2
  !  P6 (x,y,z) = x*y

  ! Corner vertice and edge midpoint coordinates
  real(DP), dimension(NDIM2D, 4) :: Dvert

  ! Coefficients for inverse affine transformation
  real(DP), dimension(NDIM2D,NDIM2D) :: Ds
  real(DP), dimension(NDIM2D) :: Dr
  real(DP) :: ddets

  ! other local variables
  integer :: iel,ipt
  real(DP) :: dx,dy

    ! Loop over all elements
    !$omp parallel do default(shared)&
    !$omp private(Dr,Ds,Dvert,ddets,dx,dy,ipt)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 1: Fetch vertice coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Fetch the four corner vertices for that element
      Dvert(1:2,1:4) = reval%p_Dcoords(1:2,1:4,iel)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 2: Calculate inverse affine transformation
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! This is a P1 transformation from the real element onto our
      ! 'reference' element.
      dr(1) = 0.25_DP * (Dvert(1,1) + Dvert(1,2) + Dvert(1,3) + Dvert(1,4))
      dr(2) = 0.25_DP * (Dvert(2,1) + Dvert(2,2) + Dvert(2,3) + Dvert(2,4))
      ds(1,1) =   0.5_DP * (Dvert(2,3) + Dvert(2,4)) - dr(2)
      ds(1,2) = -(0.5_DP * (Dvert(1,3) + Dvert(1,4)) - dr(1))
      ds(2,1) = -(0.5_DP * (Dvert(2,2) + Dvert(2,3)) - dr(2))
      ds(2,2) =   0.5_DP * (Dvert(1,2) + Dvert(1,3)) - dr(1)
      ddets = 1.0_DP / (ds(1,1)*ds(2,2) - ds(1,2)*ds(2,1))
      Ds = ddets * Ds

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 3: Evaluate function values
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_FUNC2D)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,ipt,iel) = 1.0_DP
          Dbas(2,DER_FUNC2D,ipt,iel) = dx
          DBas(3,DER_FUNC2D,ipt,iel) = dy
          Dbas(4,DER_FUNC2D,ipt,iel) = dx*dx
          DBas(5,DER_FUNC2D,ipt,iel) = dy*dy
          Dbas(6,DER_FUNC2D,ipt,iel) = dx*dy

        end do ! ipt

      end if ! function values evaluation

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Step 4: Evaluate derivatives
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

        ! Loop over all points then
        do ipt = 1, reval%npointsPerElement

          ! Apply inverse affine trafo to get (x,y)
          dx = ds(1,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(1,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))
          dy = ds(2,1)*(reval%p_DpointsReal(1,ipt,iel)-dr(1)) &
             + ds(2,2)*(reval%p_DpointsReal(2,ipt,iel)-dr(2))

          ! Evaluate derivatives
          Dbas(1,DER_DERIV2D_X,ipt,iel) = 0.0_DP
          Dbas(1,DER_DERIV2D_Y,ipt,iel) = 0.0_DP
          Dbas(2,DER_DERIV2D_X,ipt,iel) = ds(1,1)
          Dbas(2,DER_DERIV2D_Y,ipt,iel) = ds(1,2)
          Dbas(3,DER_DERIV2D_X,ipt,iel) = ds(2,1)
          Dbas(3,DER_DERIV2D_Y,ipt,iel) = ds(2,2)
          Dbas(4,DER_DERIV2D_X,ipt,iel) = 2.0_DP * dx * ds(1,1)
          Dbas(4,DER_DERIV2D_Y,ipt,iel) = 2.0_DP * dx * ds(1,2)
          Dbas(5,DER_DERIV2D_X,ipt,iel) = 2.0_DP * dy * ds(2,1)
          Dbas(5,DER_DERIV2D_Y,ipt,iel) = 2.0_DP * dy * ds(2,2)
          Dbas(6,DER_DERIV2D_X,ipt,iel) = ds(1,1)*dy + ds(2,1)*dx
          Dbas(6,DER_DERIV2D_Y,ipt,iel) = ds(1,2)*dy + ds(2,2)*dx

        end do ! ipt

      end if ! derivatives evaluation

    end do ! iel
    !$omp end parallel do

    ! That is it

  end subroutine

  !************************************************************************

!<subroutine>

 subroutine elem_eval_QPW4P1TVDF_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(in)                       :: celement

  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in)             :: reval

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in)              :: Bder
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
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! This element is piecewisely defined by 4 triangles in one
  ! quad. There are 24 DOFs: 3 DOFs on every edge of each triangle.
  !
  ! -1,1              1,1       a3                a2
  !   +-------------+             +------3------+
  !   | \    T3   / |             | \         / |
  !   |   \     /   |             |   3     3   |
  !   |     \ /     |             |     \ /     |
  !   | T4   X   T2 |     ->      3      X      3
  !   |     / \     |             |     / \     |
  !   |   /     \   |             |   3     3   |
  !   | /    T1   \ |             | /         \ |
  !   +-------------+             +------3------+
  ! -1,-1             1,-1      a0                a1
  !
  ! Therefore, for every point we have to decide in which
  ! sub-triangle we are.
  !
  ! The point ia in one of the subtriangles T1, T2, T3 or T4.
  ! Calculate the line a0->a2 and chech whether it is left
  ! or right. Calculate the line a1->a3 and check whether it is
  ! left or right. That way, we know in which subtriangle the
  ! point is.

  ! Local variables
  real(DP) :: dx,dy
  real(DP), dimension(9,9,4) :: Da
  integer :: i,iel

    ! Loop through all elements
    do iel = 1, reval%nelements
    
      ! Precalculate coefficients of the basis functions
      ! on this element.
      call elem_QPW4P1TVDF_preparetria (reval%p_Dcoords(:,:,iel),reval%p_ItwistIndex(iel),&
          1,Da(:,:,1))
      call elem_QPW4P1TVDF_preparetria (reval%p_Dcoords(:,:,iel),reval%p_ItwistIndex(iel),&
          2,Da(:,:,2))
      call elem_QPW4P1TVDF_preparetria (reval%p_Dcoords(:,:,iel),reval%p_ItwistIndex(iel),&
          3,Da(:,:,3))
      call elem_QPW4P1TVDF_preparetria (reval%p_Dcoords(:,:,iel),reval%p_ItwistIndex(iel),&
          4,Da(:,:,4))

      ! Loop through all points on the current element
      do i = 1, reval%npointsPerElement

        ! Get the point coordinates
        dx = reval%p_DpointsRef(1,i,iel)
        dy = reval%p_DpointsRef(2,i,iel)

        ! Figure out in which triangle we are.
        if (dy .le. dx) then
          ! We are in T1 or T2.
          if (dy .le. -dx) then
            ! We are in T1.
            call elem_QPW4P1TVDF_calcFunc (reval%p_Dcoords(:,:,iel),&
                reval%p_DpointsReal(:,i,iel),1,Da(:,:,1),Dbas(:,:,i,iel))
          else
            ! We are in T2
            call elem_QPW4P1TVDF_calcFunc (reval%p_Dcoords(:,:,iel),&
                reval%p_DpointsReal(:,i,iel),2,Da(:,:,2),Dbas(:,:,i,iel))
          end if
        else
          ! We are in T3 or T4
          if (dy .gt. -dx) then
            ! We are in T3
            call elem_QPW4P1TVDF_calcFunc (reval%p_Dcoords(:,:,iel),&
                reval%p_DpointsReal(:,i,iel),3,Da(:,:,3),Dbas(:,:,i,iel))
          else
            ! We are in T4
            call elem_QPW4P1TVDF_calcFunc (reval%p_Dcoords(:,:,iel),&
                reval%p_DpointsReal(:,i,iel),4,Da(:,:,4),Dbas(:,:,i,iel))
          end if
        end if

      end do ! i

    end do ! j

  end subroutine

  !************************************************************************

  subroutine elem_QPW4P1TVDF_preparetria (Dcoords,itwist,ilocaltri,Da)
  
  ! Prepares the coefficients of the basis functions on the local
  ! triangle ilocaltria 
  
  ! Array with the coordinates of the element corners.
  real(DP), dimension(:,:) :: Dcoords
  
  ! Twist index definition of the element.
  integer, intent(in) :: itwist
  
  ! Number of the local triangle.
  integer, intent(in) :: ilocaltri
  
  ! Returns the coefficients defining the local basis functions.
  real(DP), dimension(9,9), intent(out) :: Da

  ! Element Description
  ! -------------------
  ! On every sub-triangle of the quad cell, we have a set of 9 local basis 
  ! functions. Each basis function is vector valued and has the form
  !
  !     v  =  sum_{k=1,...,9}  c_k  p_k
  !
  ! where c_k is the local DOF and p_k a vector-valued polynomial defined
  ! as follows:
  !
  !    p_1 = (1,0),         p_2 = (x,0),         p_3 = (y,0)
  !    p_4 = (0,1),         p_5 = (0,x),         p_6 = (0,y)
  !    p_7 = curl(phi_1),   p_8 = curl(phi_2),   p_9 = curl(phi_3)
  !
  ! Here, the scalar-valued functions phi_i are defined by the bubble function
  ! b_T which is quadratic and =0 on the edges of the triangle:
  !
  !    phi_1 = 1 b_T,       phi_2 = x b_T,       phi_3 = y b_T
  !
  ! E.g., for the reference triangle, there is
  !
  !    b_T = x y (1-x-y)
  !
  ! But for a non-reference triangle, this is different. Assume a mapping
  !
  !    sigma: T_ref -> T,
  !
  !    sigma(xi, eta) = ( x1 ) + ( x2-x1  x3-x1 ) ( xi  )
  !                     ( y1 )   ( y2-y1  y3-y1 ) ( eta )
  !
  ! The corresponding node functionals of the element are defines using
  ! an arbitrary function v as
  !
  !    N_1(v) = 1/|E1| int_E1 (v n_1) ds
  !    N_2(v) = 1/|E2| int_E2 (v n_2) ds
  !    N_3(v) = 1/|E3| int_E3 (v n_3) ds
  !
  !    N_4(v) = 1/|E1| int_E1 (v n_1 s) ds
  !    N_5(v) = 1/|E2| int_E2 (v n_2 s) ds
  !    N_6(v) = 1/|E3| int_E3 (v n_3 s) ds
  !
  !    N_7(v) = 1/|E1| int_E1 (v t_1) ds
  !    N_8(v) = 1/|E2| int_E2 (v t_2) ds
  !    N_9(v) = 1/|E3| int_E3 (v t_3) ds
  !
  ! with n_i the normal vector of edge i and t_i the tangential.
  !
  !      x3
  !      |\
  !      | \  _ n2
  !      |  \ /|
  ! n3 <-E3 E2
  !      |    \
  !      |     \
  !      |      \
  !     x1--E1---x2
  !          |
  !          v n1
  !
  ! So,
  ! - local DOF 1..3 resemble the integral mean value of the flow through
  !   edge E1,...,E3
  ! - local DOF 4..6 resemble the 1st moment of the flow through E1,...,E3
  ! - local DOF 7..9 resemble the integral mean value of the flow along the edges.
  !
  ! The local DOFs in the triangle are ordered as follows:
  !
  !         O                               
  !         | \
  !         |  \
  !         |   \  
  !         |    \                               |          \   /           |
  !         |     \                              |            X             |
  !         |      \                             |          /   \           |
  !  3/6/9  E3~      E2~ 2/5/8                   |        /       \         |
  !         |        \              sigma        |   /  /           \  \    |
  !         |         \          ---------->     |  v / E3         E2 \ v   |
  !         |   T~     \                         |  /         T         \   |
  !         |           \                        |/                       \ |
  !         |            \                       X--------- <-o-> ----------X
  !         O-----E1~-----O                                  E1
  !             1/4/7
  !
  ! with three-tuples i/j/k representing the local DOF numbers on 
  ! edge E1,E2,E3, with
  !
  !  - i=DOF of integral mean value in normal direction
  !  - j=DOF of first moment of integral mean value in normal direction
  !  - k=DOF of integral mean value in tangential direction
  !
  ! Edge E1 corresponds to the outer edge of the triangle.
  ! Edges E2 and E3 correspond to the inner edge.
  ! By convention, the edges in the quad are oriented from the centre 
  ! to the outer, so E2 is oriented clockwise while E3 is
  ! oriented counterclockwise. The orientation of edge E1
  ! depends on the orientation of the edge in the quad and thus
  ! depends on itwistIndex.
  !
  ! The local basis functions v_1,...,v_9 corresponding to these DOFs are defined
  ! by the usual relation
  !
  !    N_i (v_j) = delta_ij
  !
  ! Thus, we have to set up a linear system for the c_k and solve it
  ! (by inverting the matrix).
  !
  ! We apply a locally nonparametric approach here, this is more easy than
  ! a parametric approach...
  
  ! Length of the edges.
  real(DP), dimension(3) :: Dlen
  
  ! X/Y coordinates of the corners of the triangle.
  real(DP), dimension(3) :: Dx, Dy
  
  ! X/Y part of the tangential and normal vectors on the edges.
  real(DP), dimension(3) :: Dtx, Dty, Dnx, Dny
  
  ! Parameter values and weights of cubature points on each edge.
  ! We apply a 3-point Gauss formula on every edge for setting
  ! up the N_i(p_j).
  integer, parameter :: ncubp = 3
  
  real(DP), dimension(ncubp), parameter :: Dg = &
    (/ -0.774596669241483_DP, 0.0_DP, 0.774596669241483_DP /)

  real(DP), dimension(ncubp), parameter :: Dw = &
    (/ 0.555555555555556_DP, 0.888888888888889_DP, 0.555555555555556_DP /)
    
  ! X/Y position of the cubature points on each edge
  real(DP), dimension(ncubp,3) :: DcubpXref,DcubpYref
  
  ! Curl of the trial polynomial
  real(DP),dimension(ncubp,3) :: Dcurlx1,Dcurly1
  real(DP),dimension(ncubp,3) :: Dcurlx2,Dcurly2
  real(DP),dimension(ncubp,3) :: Dcurlx3,Dcurly3
  real(DP), dimension(ncubp) :: Dlambda
  
  real(DP) :: dtwist1,dtwist2,dtwist3
  real(DP) :: ddet,a11,a12,a21,a22
  integer :: i,j
  logical :: bsuccess

    ! Get the twist indices that define the orientation of the four edges
    ! A value of 1 is standard, a value of -1 results in a change of the
    ! sign in the basis functions.
    !
    ! Edge 1 corresponds to the outer edge of the quad.
    ! Edge 2/3 corresponds to the inner edges. By convention, the
    ! edges in the quad are oriented from the centre to the outer, so
    ! the twist index of the 2nd edge is always -1 and the twist index
    ! of the 3rd edge always 1.
    !
    ! itwistIndex is a bitfield. Each bit specifies the orientation of an edge.
    ! We use the bit to calculate a "1.0" if the edge has positive orientation
    ! and "-1.0" if it has negative orientation.
    dtwist1 = real(1-iand(int(ishft(itwist, 2-ilocaltri)),2),DP)
    dtwist2 = -1.0_DP
    dtwist3 = 1.0_DP

    ! Local Element corners
    Dx(1) = Dcoords(1,ilocaltri)
    Dy(1) = Dcoords(2,ilocaltri)
    Dx(2) = Dcoords(1,mod(ilocaltri,4)+1)
    Dy(2) = Dcoords(2,mod(ilocaltri,4)+1)
    Dx(3) = 0.25_DP*sum(Dcoords(1,1:4))
    Dy(3) = 0.25_DP*sum(Dcoords(2,1:4))

    ! Transformation from the reference to the real element
    a11 = Dx(2)-Dx(1)
    a12 = Dx(3)-Dx(1)
    a21 = Dy(2)-Dy(1)
    a22 = Dy(3)-Dy(1)
    ddet = 1.0_DP / (a11*a22-a12*a21)

    ! Edge lengths
    Dlen(1) = sqrt( (Dx(2)-Dx(1))**2 + (Dy(2)-Dy(1))**2 )
    Dlen(2) = sqrt( (Dx(3)-Dx(2))**2 + (Dy(3)-Dy(2))**2 )
    Dlen(3) = sqrt( (Dx(1)-Dx(3))**2 + (Dy(1)-Dy(3))**2 )
    
    ! Calculate the tangential and normal vectors to the edges from the corners.
    ! The tangential vectors of the edges are
    !
    !    t1 = (x2-x1) / ||x2-x1||
    !    t2 = (x3-x2) / ||x3-x2||
    !    t3 = (x1-x3) / ||x1-x3||
    !
    ! The corresponding normals are calculated from
    !
    !    ni = ( ti2, -ti1 )
    !
    Dtx(1) = (Dx(2)-Dx(1)) / Dlen(1) * dtwist1
    Dty(1) = (Dy(2)-Dy(1)) / Dlen(1) * dtwist1

    Dtx(2) = (Dx(3)-Dx(2)) / Dlen(2) * dtwist2
    Dty(2) = (Dy(3)-Dy(2)) / Dlen(2) * dtwist2

    Dtx(3) = (Dx(1)-Dx(3)) / Dlen(3) * dtwist3
    Dty(3) = (Dy(1)-Dy(3)) / Dlen(3) * dtwist3

    Dnx(1) = Dty(1)
    Dny(1) = -Dtx(1)

    Dnx(2) = Dty(2)
    Dny(2) = -Dtx(2)

    Dnx(3) = Dty(3)
    Dny(3) = -Dtx(3)
    
    ! Set up a transformation matrix to calculate the c_k.
    !
    ! The conditions to be solved for one basis function read:
    !
    !   N_i(phi_j) = delta_ij
    !
    ! For the calculation of the N_i, we apply a point Gauss
    ! formula. The basic coordinates of the Gauss formula on
    ! [-1,1] and the cubature weights are defined as
    
    ! Calculate the cubature points for the evaluation of the coefficients.
    ! This is on the reference element.
    !
    ! Orientation of the first edge must possibly be inverted,
    ! depending on the twist index.
    if (dtwist1 .gt. 0.0_DP) then
      do i=1,ncubp
        ! Parameter value in [0,1]
        Dlambda(i) = 0.5_DP * (Dg(i) + 1.0_DP)

        ! Cubature points on the reference triangle.
        DcubpXref(i,1) = Dlambda(i) ! * 1.0 + (1.0_DP - Dlambda(i)) * 0.0_DP
        DcubpYref(i,1) = 0.0_DP
        
        ! Orientation of the second edge must be inverted.
        DcubpXref(i,2) = Dlambda(i) ! * 1.0_DP + (1.0_DP - Dlambda(i)) * 0.0_DP 
        DcubpYref(i,2) = (1.0_DP - Dlambda(i)) ! * 1.0_DP + Dlambda(i) * 0.0_DP + 
        
        DcubpXref(i,3) = 0.0_DP
        DcubpYref(i,3) = (1.0_DP - Dlambda(i)) ! * 1.0_DP + Dlambda(i) * 0.0
      end do
    else
      do i=1,ncubp
        ! Parameter value in [0,1]
        Dlambda(i) = 0.5_DP * (Dg(i) + 1.0_DP)

        ! Cubature points on the reference triangle.
        DcubpXref(i,1) = 1.0_DP - Dlambda(i)
        DcubpYref(i,1) = 0.0_DP
        
        ! Orientation of the second edge must be inverted.
        DcubpXref(i,2) = Dlambda(i) ! * 1.0_DP + (1.0_DP - Dlambda(i)) * 0.0_DP 
        DcubpYref(i,2) = (1.0_DP - Dlambda(i)) ! * 1.0_DP + Dlambda(i) * 0.0_DP + 
        
        DcubpXref(i,3) = 0.0_DP
        DcubpYref(i,3) = (1.0_DP - Dlambda(i)) ! * 1.0_DP + Dlambda(i) * 0.0
      end do
    end if

    ! Calculate the matrix N=(n_ij) with
    !
    !    n_ij = N_i(p_j)
    !
    ! Note that the normal vector is constant along the complete edge,
    ! which simplifies this a bit. Furthermore note that
    ! the length of the edge (in the cubature weight) cancels out
    ! against the normalisation factor in front of the integral,
    ! which simplifies this even more.
    
    ! N_i(p_1)
    Da(1,1) = 0.5_DP * sum ( Dnx(1) * Dw(:) )
    Da(2,1) = 0.5_DP * sum ( Dnx(2) * Dw(:) )
    Da(3,1) = 0.5_DP * sum ( Dnx(3) * Dw(:) )
    Da(4,1) = 0.5_DP * sum ( Dnx(1) * Dg(:) * Dw(:) )
    Da(5,1) = 0.5_DP * sum ( Dnx(2) * Dg(:) * Dw(:) )
    Da(6,1) = 0.5_DP * sum ( Dnx(3) * Dg(:) * Dw(:) )
    Da(7,1) = 0.5_DP * sum ( Dtx(1) * Dw(:) )
    Da(8,1) = 0.5_DP * sum ( Dtx(2) * Dw(:) )
    Da(9,1) = 0.5_DP * sum ( Dtx(3) * Dw(:) )
    
    ! N_i(p_2)
    Da(1,2) = 0.5_DP * sum ( DcubpXref(:,1) * Dnx(1) * Dw(:) )
    Da(2,2) = 0.5_DP * sum ( DcubpXref(:,2) * Dnx(2) * Dw(:) )
    Da(3,2) = 0.5_DP * sum ( DcubpXref(:,3) * Dnx(3) * Dw(:) )
    Da(4,2) = 0.5_DP * sum ( DcubpXref(:,1) * Dnx(1) * Dg(:) * Dw(:) )
    Da(5,2) = 0.5_DP * sum ( DcubpXref(:,2) * Dnx(2) * Dg(:) * Dw(:) )
    Da(6,2) = 0.5_DP * sum ( DcubpXref(:,3) * Dnx(3) * Dg(:) * Dw(:) )
    Da(7,2) = 0.5_DP * sum ( DcubpXref(:,1) * Dtx(1) * Dw(:) )
    Da(8,2) = 0.5_DP * sum ( DcubpXref(:,2) * Dtx(2) * Dw(:) )
    Da(9,2) = 0.5_DP * sum ( DcubpXref(:,3) * Dtx(3) * Dw(:) )

    ! N_i(p_3)
    Da(1,3) = 0.5_DP * sum ( DcubpYref(:,1) * Dnx(1) * Dw(:) )
    Da(2,3) = 0.5_DP * sum ( DcubpYref(:,2) * Dnx(2) * Dw(:) )
    Da(3,3) = 0.5_DP * sum ( DcubpYref(:,3) * Dnx(3) * Dw(:) )
    Da(4,3) = 0.5_DP * sum ( DcubpYref(:,1) * Dnx(1) * Dg(:) * Dw(:) )
    Da(5,3) = 0.5_DP * sum ( DcubpYref(:,2) * Dnx(2) * Dg(:) * Dw(:) )
    Da(6,3) = 0.5_DP * sum ( DcubpYref(:,3) * Dnx(3) * Dg(:) * Dw(:) )
    Da(7,3) = 0.5_DP * sum ( DcubpYref(:,1) * Dtx(1) * Dw(:) )
    Da(8,3) = 0.5_DP * sum ( DcubpYref(:,2) * Dtx(2) * Dw(:) )
    Da(9,3) = 0.5_DP * sum ( DcubpYref(:,3) * Dtx(3) * Dw(:) )
    
    ! N_i(p_4)
    Da(1,4) = 0.5_DP * sum ( Dny(1) * Dw(:) )
    Da(2,4) = 0.5_DP * sum ( Dny(2) * Dw(:) )
    Da(3,4) = 0.5_DP * sum ( Dny(3) * Dw(:) )
    Da(4,4) = 0.5_DP * sum ( Dny(1) * Dg(:) * Dw(:) )
    Da(5,4) = 0.5_DP * sum ( Dny(2) * Dg(:) * Dw(:) )
    Da(6,4) = 0.5_DP * sum ( Dny(3) * Dg(:) * Dw(:) )
    Da(7,4) = 0.5_DP * sum ( Dty(1) * Dw(:) )
    Da(8,4) = 0.5_DP * sum ( Dty(2) * Dw(:) )
    Da(9,4) = 0.5_DP * sum ( Dty(3) * Dw(:) )
    
    ! N_i(p_5)
    Da(1,5) = 0.5_DP * sum ( DcubpXref(:,1) * Dny(1) * Dw(:) )
    Da(2,5) = 0.5_DP * sum ( DcubpXref(:,2) * Dny(2) * Dw(:) )
    Da(3,5) = 0.5_DP * sum ( DcubpXref(:,3) * Dny(3) * Dw(:) )
    Da(4,5) = 0.5_DP * sum ( DcubpXref(:,1) * Dny(1) * Dg(:) * Dw(:) )
    Da(5,5) = 0.5_DP * sum ( DcubpXref(:,2) * Dny(2) * Dg(:) * Dw(:) )
    Da(6,5) = 0.5_DP * sum ( DcubpXref(:,3) * Dny(3) * Dg(:) * Dw(:) )
    Da(7,5) = 0.5_DP * sum ( DcubpXref(:,1) * Dty(1) * Dw(:) )
    Da(8,5) = 0.5_DP * sum ( DcubpXref(:,2) * Dty(2) * Dw(:) )
    Da(9,5) = 0.5_DP * sum ( DcubpXref(:,3) * Dty(3) * Dw(:) )
    
    ! N_i(p_6)
    Da(1,6) = 0.5_DP * sum ( DcubpYref(:,1) * Dny(1) * Dw(:) )
    Da(2,6) = 0.5_DP * sum ( DcubpYref(:,2) * Dny(2) * Dw(:) )
    Da(3,6) = 0.5_DP * sum ( DcubpYref(:,3) * Dny(3) * Dw(:) )
    Da(4,6) = 0.5_DP * sum ( DcubpYref(:,1) * Dny(1) * Dg(:) * Dw(:) )
    Da(5,6) = 0.5_DP * sum ( DcubpYref(:,2) * Dny(2) * Dg(:) * Dw(:) )
    Da(6,6) = 0.5_DP * sum ( DcubpYref(:,3) * Dny(3) * Dg(:) * Dw(:) )
    Da(7,6) = 0.5_DP * sum ( DcubpYref(:,1) * Dty(1) * Dw(:) )
    Da(8,6) = 0.5_DP * sum ( DcubpYref(:,2) * Dty(2) * Dw(:) )
    Da(9,6) = 0.5_DP * sum ( DcubpYref(:,3) * Dty(3) * Dw(:) )

    ! Evaluate P7, P8, P9 in the cubature points.
    do j=1,3
      do i=1,ncubp
        call QPW4P1TVDF_P1(DcubpXref(i,j), DcubpYref(i,j), &
            a11,a12,a21,a22,ddet, Dcurlx1(i,j),Dcurly1(i,j))
        call QPW4P1TVDF_P2(DcubpXref(i,j), DcubpYref(i,j), &
            a11,a12,a21,a22,ddet, Dcurlx2(i,j),Dcurly2(i,j))
        call QPW4P1TVDF_P3(DcubpXref(i,j), DcubpYref(i,j), &
            a11,a12,a21,a22,ddet, Dcurlx3(i,j),Dcurly3(i,j))
      end do
    end do

    ! N_i(p_7)
    Da(1,7) = 0.5_DP * sum ( Dnx(1) * Dcurlx1(:,1) * Dw(:) &
                           + Dny(1) * Dcurly1(:,1) * Dw(:) )
    Da(2,7) = 0.5_DP * sum ( Dnx(2) * Dcurlx1(:,2) * Dw(:) &
                           + Dny(2) * Dcurly1(:,2) * Dw(:) )
    Da(3,7) = 0.5_DP * sum ( Dnx(3) * Dcurlx1(:,3) * Dw(:) &
                           + Dny(3) * Dcurly1(:,3) * Dw(:) )
    Da(4,7) = 0.5_DP * sum ( Dnx(1) * Dcurlx1(:,1) * Dg(:) * Dw(:) &
                           + Dny(1) * Dcurly1(:,1) * Dg(:) * Dw(:) )
    Da(5,7) = 0.5_DP * sum ( Dnx(2) * Dcurlx1(:,2) * Dg(:) * Dw(:) &
                           + Dny(2) * Dcurly1(:,2) * Dg(:) * Dw(:) )
    Da(6,7) = 0.5_DP * sum ( Dnx(3) * Dcurlx1(:,3) * Dg(:) * Dw(:) &
                           + Dny(3) * Dcurly1(:,3) * Dg(:) * Dw(:) )
    Da(7,7) = 0.5_DP * sum ( Dtx(1) * Dcurlx1(:,1) * Dw(:) &
                           + Dty(1) * Dcurly1(:,1) * Dw(:) )
    Da(8,7) = 0.5_DP * sum ( Dtx(2) * Dcurlx1(:,2) * Dw(:) &
                           + Dty(2) * Dcurly1(:,2) * Dw(:) )
    Da(9,7) = 0.5_DP * sum ( Dtx(3) * Dcurlx1(:,3) * Dw(:) &
                           + Dty(3) * Dcurly1(:,3) * Dw(:) )

    ! N_i(p_8)
    Da(1,8) = 0.5_DP * sum ( Dnx(1) * Dcurlx2(:,1) * Dw(:) &
                           + Dny(1) * Dcurly2(:,1) * Dw(:) )
    Da(2,8) = 0.5_DP * sum ( Dnx(2) * Dcurlx2(:,2) * Dw(:) &
                           + Dny(2) * Dcurly2(:,2) * Dw(:) )
    Da(3,8) = 0.5_DP * sum ( Dnx(3) * Dcurlx2(:,3) * Dw(:) &
                           + Dny(3) * Dcurly2(:,3) * Dw(:) )
    Da(4,8) = 0.5_DP * sum ( Dnx(1) * Dcurlx2(:,1) * Dg(:) * Dw(:) &
                           + Dny(1) * Dcurly2(:,1) * Dg(:) * Dw(:) )
    Da(5,8) = 0.5_DP * sum ( Dnx(2) * Dcurlx2(:,2) * Dg(:) * Dw(:) &
                           + Dny(2) * Dcurly2(:,2) * Dg(:) * Dw(:) )
    Da(6,8) = 0.5_DP * sum ( Dnx(3) * Dcurlx2(:,3) * Dg(:) * Dw(:) &
                           + Dny(3) * Dcurly2(:,3) * Dg(:) * Dw(:) )
    Da(7,8) = 0.5_DP * sum ( Dtx(1) * Dcurlx2(:,1) * Dw(:) &
                           + Dty(1) * Dcurly2(:,1) * Dw(:) )
    Da(8,8) = 0.5_DP * sum ( Dtx(2) * Dcurlx2(:,2) * Dw(:) &
                           + Dty(2) * Dcurly2(:,2) * Dw(:) )
    Da(9,8) = 0.5_DP * sum ( Dtx(3) * Dcurlx2(:,3) * Dw(:) &
                           + Dty(3) * Dcurly2(:,3) * Dw(:) )
    
    ! N_i(p_9)
    Da(1,9) = 0.5_DP * sum ( Dnx(1) * Dcurlx3(:,1) * Dw(:) &
                           + Dny(1) * Dcurly3(:,1) * Dw(:) )
    Da(2,9) = 0.5_DP * sum ( Dnx(2) * Dcurlx3(:,2) * Dw(:) &
                           + Dny(2) * Dcurly3(:,2) * Dw(:) )
    Da(3,9) = 0.5_DP * sum ( Dnx(3) * Dcurlx3(:,3) * Dw(:) &
                           + Dny(3) * Dcurly3(:,3) * Dw(:) )
    Da(4,9) = 0.5_DP * sum ( Dnx(1) * Dcurlx3(:,1) * Dg(:) * Dw(:) &
                           + Dny(1) * Dcurly3(:,1) * Dg(:) * Dw(:) )
    Da(5,9) = 0.5_DP * sum ( Dnx(2) * Dcurlx3(:,2) * Dg(:) * Dw(:) &
                           + Dny(2) * Dcurly3(:,2) * Dg(:) * Dw(:) )
    Da(6,9) = 0.5_DP * sum ( Dnx(3) * Dcurlx3(:,3) * Dg(:) * Dw(:) &
                           + Dny(3) * Dcurly3(:,3) * Dg(:) * Dw(:) )
    Da(7,9) = 0.5_DP * sum ( Dtx(1) * Dcurlx3(:,1) * Dw(:) &
                           + Dty(1) * Dcurly3(:,1) * Dw(:) )
    Da(8,9) = 0.5_DP * sum ( Dtx(2) * Dcurlx3(:,2) * Dw(:) &
                           + Dty(2) * Dcurly3(:,2) * Dw(:) )
    Da(9,9) = 0.5_DP * sum ( Dtx(3) * Dcurlx3(:,3) * Dw(:) &
                           + Dty(3) * Dcurly3(:,3) * Dw(:) )

    ! Taking the inverse will give us the coefficients of all
    ! basis functions:
    !
    !    A C = Id  =>  C=Inverse of the matrix
    call mprim_invertMatrixPivot(Da,9,bsuccess)
    
  end subroutine
    
  !************************************************************************

  subroutine elem_QPW4P1TVDF_calcFunc (Dcoords,Dpoint,ilocaltri,Da,Dbas)
  
  ! Calculates the values of the local basis functions in point Ipoint. 
  
  ! Corners of the element
  real(DP), dimension(:,:), intent(in) :: Dcoords
  
  ! Coordinate of the point where to evaluate
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Number of the local triangle
  integer, intent(in) :: ilocaltri
  
  ! Coefficients defining the local basis functions on triangle ilocaltria
  real(DP), dimension(:,:), intent(in) :: Da
  
  ! Returns the values of the local basis functions in Dpoint.
  real(DP), dimension(:,:), intent(out) :: Dbas
  
  ! Local variables
  real(DP), dimension(3) :: Dx,Dy
  real(DP), dimension(9) :: DbasisX, DbasisY
  real(DP), dimension(9) :: DbasisX_X, DbasisY_X
  real(DP), dimension(9) :: DbasisX_Y, DbasisY_Y
  real(DP) :: dpx, dpy, ddet
  real(DP) :: dx_x,dx_y,dy_x,dy_y, a11,a12,a21,a22
  integer :: i
  
  ! The following array defines the mapping of the DOFs from the
  ! local triangle to the local quad. The DOFs in the quad are
  ! ordered as follows:
  !
  !                        3/7/11
  !      O-------------------X-------------------O
  !      | \                                   / |
  !      |   \                               /   |
  !      |     \                           /     |
  !      |       \  16/20/24             /       |
  !      |         X                   X         |
  !      |           \               /  15/19/23 |
  !      |             \           /             |
  !      |               \       /               |
  !      |                 \   /                 |
  !      X 4/8/12            O                   X 2/6/10
  !      |                 /   \                 |
  !      |               /       \               |
  !      |             /           \             |
  !      |           /               \  14/18/22 |
  !      |         X                   X         |
  !      |       /  13/17/21             \       |
  !      |     /                           \     |
  !      |   /                               \   |
  !      | /                                   \ |
  !      O-------------------X-------------------O
  !                        1/5/9
  !
  ! Again, the three-tuples i/j/k representing the local DOF
  !
  !  - i=DOF of integral mean value in normal direction
  !  - j=DOF of first moment of integral mean value in normal direction
  !  - k=DOF of integral mean value in tangential direction
  !
  integer, dimension(9,4), parameter :: Ildofs = reshape (&
    (/ 1, 14, 13, 5, 18, 17,  9, 22, 21, &
       2, 15, 14, 6, 19, 18, 10, 23, 22, &
       3, 16, 15, 7, 20, 19, 11, 24, 23, &
       4, 13, 16, 8, 17, 20, 12, 21, 24 /), (/ 9, 4 /) )

    ! Clear the output array.
    Dbas(:,DER_FUNC2D) = 0.0_DP
    Dbas(:,DER_DERIV2D_X) = 0.0_DP
    Dbas(:,DER_DERIV2D_Y) = 0.0_DP
    
    ! Local Element corners
    Dx(1) = Dcoords(1,ilocaltri)
    Dy(1) = Dcoords(2,ilocaltri)
    Dx(2) = Dcoords(1,mod(ilocaltri,4)+1)
    Dy(2) = Dcoords(2,mod(ilocaltri,4)+1)
    Dx(3) = 0.25_DP*sum(Dcoords(1,1:4))
    Dy(3) = 0.25_DP*sum(Dcoords(2,1:4))
    
    ! Transform the point to the reference triangle.
    a11 = Dx(2)-Dx(1)
    a12 = Dx(3)-Dx(1)
    a21 = Dy(2)-Dy(1)
    a22 = Dy(3)-Dy(1)
    ddet = 1.0_DP / (a11*a22 - a12*a21)
    dpx = ddet * ( ( a22) * (Dpoint(1)-Dx(1)) +  (-a12) * (Dpoint(2) - Dy(1)) )
    dpy = ddet * ( (-a21) * (Dpoint(1)-Dx(1)) +  ( a11) * (Dpoint(2) - Dy(1)) )

    ! Evaluate the in the point; function values.
    DbasisX(1) = 1.0_DP
    DbasisX(2) = dpx
    DbasisX(3) = dpy
    DbasisX(4) = 0.0_DP
    DbasisX(5) = 0.0_DP
    DbasisX(6) = 0.0_DP

    DbasisY(1) = 0.0_DP
    DbasisY(2) = 0.0_DP
    DbasisY(3) = 0.0_DP
    DbasisY(4) = 1.0_DP
    DbasisY(5) = dpx
    DbasisY(6) = dpy

    call QPW4P1TVDF_P1(dpx, dpy, a11,a12,a21,a22,ddet, DbasisX(7), DbasisY(7))
    call QPW4P1TVDF_P2(dpx, dpy, a11,a12,a21,a22,ddet, DbasisX(8), DbasisY(8))
    call QPW4P1TVDF_P3(dpx, dpy, a11,a12,a21,a22,ddet, DbasisX(9), DbasisY(9))

    ! Evaluate basis functions:
    !
    !   v = sum_{i=1,9} x_i p_i
    
    do i=1,9
      ! Dofs 1..24 are the 1st coordinate, DOFs 25..48 the 2nd coordinate
      ! of the vector-valued element.

      ! Function value.
      Dbas (Ildofs(i,ilocaltri)   ,DER_FUNC2D) = sum ( Da(:,i) * DbasisX(:) )
      Dbas (Ildofs(i,ilocaltri)+24,DER_FUNC2D) = sum ( Da(:,i) * DbasisY(:) )

    end do
            
    ! Evaluate the in the point; derivatives
    DbasisX_X(1) = 0.0_DP
    DbasisX_X(2) = 1.0_DP
    DbasisX_X(3) = 0.0_DP
    DbasisX_X(4) = 0.0_DP
    DbasisX_X(5) = 0.0_DP
    DbasisX_X(6) = 0.0_DP

    DbasisY_X(1) = 0.0_DP
    DbasisY_X(2) = 0.0_DP
    DbasisY_X(3) = 0.0_DP
    DbasisY_X(4) = 0.0_DP
    DbasisY_X(5) = 1.0_DP
    DbasisY_X(6) = 0.0_DP

    DbasisX_Y(1) = 0.0_DP
    DbasisX_Y(2) = 0.0_DP
    DbasisX_Y(3) = 1.0_DP
    DbasisX_Y(4) = 0.0_DP
    DbasisX_Y(5) = 0.0_DP
    DbasisX_Y(6) = 0.0_DP

    DbasisY_Y(1) = 0.0_DP
    DbasisY_Y(2) = 0.0_DP
    DbasisY_Y(3) = 0.0_DP
    DbasisY_Y(4) = 0.0_DP
    DbasisY_Y(5) = 0.0_DP
    DbasisY_Y(6) = 1.0_DP

    call QPW4P1TVDF_P1_XY(dpx, dpy, &
        a11,a12,a21,a22,ddet, DbasisX_X(7), DbasisY_X(7), DbasisX_Y(7), DbasisY_Y(7))
    call QPW4P1TVDF_P2_XY(dpx, dpy, &
        a11,a12,a21,a22,ddet, DbasisX_X(8), DbasisY_X(8), DbasisX_Y(8), DbasisY_Y(8))
    call QPW4P1TVDF_P3_XY(dpx, dpy, &
        a11,a12,a21,a22,ddet, DbasisX_X(9), DbasisY_X(9), DbasisX_Y(9), DbasisY_Y(9))

    ! Evaluate first derivative.
    !
    !   v_x = d/dx sum_{i=1,9} x_i p_i
    !   v_y = d/dy sum_{i=1,9} x_i p_i
    
    do i=1,9
      ! Dofs 1..24 are the 1st coordinate, DOFs 25..48 the 2nd coordinate
      ! of the vector-valued element.

      ! Derivatives values on the reference element
      !
      ! Multiply by the inverse of the jacobian mapping from the
      ! reference triangle to the real triangle in order to
      ! calculate the actual derivative.
      !
      ! Every basis function is defined as
      !
      !    phi_i = sum_j c_ij p_j
      !
      ! with p_j the vector valued trial polynomials defined as
      !
      !    p_j(x) = p^_j ( sigma^-1 (x) )
      !
      ! and sigma the mapping from the reference triangle to the
      ! local triangle ilocaltri in the quad, given as
      !
      !    sigma(x^) = ( a11 a12 ) x^ + ( x1 )
      !                ( a21 a22 )      ( y1 )
      !
      ! Due to the chain rule, the derivative of phi_i is given by
      !
      !    D phi (x) = D (p^_j ( sigma^-1 (x) ) )
      !              = D p^_j (x^) * [D sigma(x)]^-1
      !
      ! The above dx_x,...,dy_y define "D p^_j (x^)". Multiply this
      ! with [D sigma(x)]^-1 to get "D phi(x)".

      dx_x = sum ( Da(:,i) * ( DbasisX_X(:)*a22 - DbasisX_Y(:)*a21 ) )
      dy_x = sum ( Da(:,i) * ( DbasisY_X(:)*a22 - DbasisY_Y(:)*a21 ) )

      dx_y = sum ( Da(:,i) * (-DbasisX_X(:)*a12 + DbasisX_Y(:)*a11 ) )
      dy_y = sum ( Da(:,i) * (-DbasisY_X(:)*a12 + DbasisY_Y(:)*a11 ) )

      Dbas (Ildofs(i,ilocaltri)   ,DER_DERIV2D_X) = ddet*dx_x
      Dbas (Ildofs(i,ilocaltri)+24,DER_DERIV2D_X) = ddet*dy_x
                                                         
      Dbas (Ildofs(i,ilocaltri)   ,DER_DERIV2D_Y) = ddet*dx_y
      Dbas (Ildofs(i,ilocaltri)+24,DER_DERIV2D_Y) = ddet*dy_y
    end do
    
  end subroutine
  
  subroutine QPW4P1TVDF_P1 (dx,dy,a11,a12,a21,a22,ddet,&
      dcurlx,dcurly)
  real(DP), intent(in) :: dx,dy
  real(DP), intent(in) :: a11,a12,a21,a22,ddet
  real(DP), intent(out) :: dcurlx, dcurly
  
  real(DP) :: x1x,x1y,x2x,x2y
    ! P1 = curl (b_T) ... plus chain rule, as this is realised on the real element
    ! dcurlx = dx-dx**2-2*dx*dy
    ! dcurly = -dy+2*dx*dy+dy**2
    x1x = a22*ddet
    x1y = -a12*ddet
    x2x = -a21*ddet
    x2y = a11*ddet

    dcurlx = x1y*dy-2*x1y*dy*dx-x1y*dy**2+dx*x2y-dx**2*x2y-2*dx*dy*x2y
    dcurly = -x1x*dy+2*dx*dy*x1x+x1x*dy**2-dx*x2x+dx**2*x2x+2*dx*dy*x2x
    
  end subroutine
  
  subroutine QPW4P1TVDF_P2 (dx,dy,a11,a12,a21,a22,ddet,&
      dcurlx,dcurly)
  real(DP), intent(in) :: dx,dy
  real(DP), intent(in) :: a11,a12,a21,a22,ddet
  real(DP), intent(out) :: dcurlx, dcurly
  
  real(DP) :: x1x,x1y,x2x,x2y
    ! P2 = curl (b_T x) ... plus chain rule, as this is realised on the real element
    ! dcurlx = dx**2-dx**3-2*dx**2*dy
    ! dcurly = -2*dx*dy+3*dx**2*dy+2*dx*dy**2

    x1x = a22*ddet
    x1y = -a12*ddet
    x2x = -a21*ddet
    x2y = a11*ddet

    dcurlx = 2*x1y*dy*dx-3*x1y*dy*dx**2-2*x1y*dy**2*dx+dx**2*x2y-dx**3*x2y-2*dx**2*x2y*dy
    dcurly = -2*dx*dy*x1x+3*x1x*dy*dx**2+2*x1x*dy**2*dx-dx**2*x2x+dx**3*x2x+2*dx**2*x2x*dy

  end subroutine

  subroutine QPW4P1TVDF_P3 (dx,dy,a11,a12,a21,a22,ddet,&
      dcurlx,dcurly)
  real(DP), intent(in) :: dx,dy
  real(DP), intent(in) :: a11,a12,a21,a22,ddet
  real(DP), intent(out) :: dcurlx, dcurly
  
  real(DP) :: x1x,x1y,x2x,x2y
    ! P3 = curl (b_T * y) ... plus chain rule, as this is realised on the real element
    ! dcurlx = 2*dx*dy-2*dx**2*dy-3*dx*dy**2
    ! dcurly = -dy**2+2*dx*dy**2+dy**3

    x1x = a22*ddet
    x1y = -a12*ddet
    x2x = -a21*ddet
    x2y = a11*ddet

    dcurlx = x1y*dy**2-2*x1y*dy**2*dx-x1y*dy**3+2*dx*dy*x2y-2*dx**2*x2y*dy-3*dx*x2y*dy**2
    dcurly = -x1x*dy**2+2*x1x*dy**2*dx+x1x*dy**3-2*dx*dy*x2x+2*dx**2*x2x*dy+3*dx*dy**2*x2x
  end subroutine

  ! 1st derivative of the curl functions
  subroutine QPW4P1TVDF_P1_XY (dx,dy,a11,a12,a21,a22,ddet,&
      dcurlx_x,dcurly_x,dcurlx_y,dcurly_y)
  real(DP), intent(in) :: dx,dy
  real(DP), intent(in) :: a11,a12,a21,a22,ddet
  real(DP), intent(out) :: dcurlx_x, dcurly_x, dcurlx_y, dcurly_y
  
  real(DP) :: x1x,x1y,x2x,x2y
    ! P1 = curl (b_T)

    x1x = a22*ddet
    x1y = -a12*ddet
    x2x = -a21*ddet
    x2y = a11*ddet

    dcurlx_x = -2*x1y*dy+x2y-2*dx*x2y-2*dy*x2y
    dcurlx_y = x1y-2*x1y*dx-2*x1y*dy-2*dx*x2y
    dcurly_x = 2*x1x*dy-x2x+2*dx*x2x+2*dy*x2x
    dcurly_y = -x1x+2*dx*x1x+2*x1x*dy+2*dx*x2x
  end subroutine
  
  subroutine QPW4P1TVDF_P2_XY (dx,dy,a11,a12,a21,a22,ddet,&
      dcurlx_x,dcurly_x,dcurlx_y,dcurly_y)
  real(DP), intent(in) :: dx,dy
  real(DP), intent(in) :: a11,a12,a21,a22,ddet
  real(DP), intent(out) :: dcurlx_x, dcurly_x, dcurlx_y, dcurly_y
  
  real(DP) :: x1x,x1y,x2x,x2y
    ! P2 = curl (b_T x)

    x1x = a22*ddet
    x1y = -a12*ddet
    x2x = -a21*ddet
    x2y = a11*ddet
    
    dcurlx_x = 2*x1y*dy-6*x1y*dy*dx-2*x1y*dy**2+2*dx*x2y-3*dx**2*x2y-4*dx*dy*x2y
    dcurlx_y = 2*x1y*dx-3*x1y*dx**2-4*x1y*dy*dx-2*dx**2*x2y
    dcurly_x = -2*x1x*dy+6*dx*dy*x1x+2*x1x*dy**2-2*dx*x2x+3*dx**2*x2x+4*dx*dy*x2x
    dcurly_y = -2*dx*x1x+3*x1x*dx**2+4*dx*dy*x1x+2*dx**2*x2x
  end subroutine

  subroutine QPW4P1TVDF_P3_XY (dx,dy,a11,a12,a21,a22,ddet,&
      dcurlx_x,dcurly_x,dcurlx_y,dcurly_y)
  real(DP), intent(in) :: dx,dy
  real(DP), intent(in) :: a11,a12,a21,a22,ddet
  real(DP), intent(out) :: dcurlx_x, dcurly_x, dcurlx_y, dcurly_y
  
  real(DP) :: x1x,x1y,x2x,x2y
    ! P3 = curl (b_T * y)

    x1x = a22*ddet
    x1y = -a12*ddet
    x2x = -a21*ddet
    x2y = a11*ddet
    
    dcurlx_x = -2*x1y*dy**2+2*dy*x2y-4*dx*dy*x2y-3*x2y*dy**2
    dcurlx_y = 2*x1y*dy-4*x1y*dy*dx-3*x1y*dy**2+2*dx*x2y-2*dx**2*x2y-6*dx*dy*x2y
    dcurly_x = 2*x1x*dy**2-2*dy*x2x+4*dx*dy*x2x+3*dy**2*x2x
    dcurly_y = -2*x1x*dy+4*dx*dy*x1x+3*x1x*dy**2-2*dx*x2x+2*dx**2*x2x+6*dx*dy*x2x
  end subroutine

end module
