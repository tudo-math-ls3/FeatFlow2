!##############################################################################
!# ****************************************************************************
!# <name> element_tri2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the implementations of the 2D triangle basis functions.
!#
!# </purpose>
!##############################################################################

module element_tri2d

!$use omp_lib
  use fsystem
  use elementbase
  use derivatives
  use perfconfig
  use mprimitives

  implicit none

  private

  public :: elem_P0
  public :: elem_P0_mult
  public :: elem_P0_sim
  public :: elem_P1
  public :: elem_P1_mult
  public :: elem_P1_sim
  public :: elem_P2
  public :: elem_P2_mult
  public :: elem_P2_sim
  public :: elem_P1T
  public :: elem_P1T_mult
  public :: elem_P1T_sim
  
  public :: elem_eval_RT1_2D
  public :: elem_eval_DCP1_2D

contains

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
! Let us use the notation
!
!    P(X,Y)  :=  p( t^-1 (X,Y) )  =  p (xi1,xi2,xi3)
!
! which is the 'standard' form of the polynomials on the
! reference element without barycentric coordinates involved (i.e.
! for the P1-element it is of the form P(x,y) = c1 x + c2 y + c3).
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

  pure subroutine elem_P0 (celement, Dcoords, Djac, ddetj, Bder, &
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

  ! Barycentric coordinates of the point where to evaluate
  real(DP), dimension(3), intent(in) :: Dpoint
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

  ! We have no derivatives.

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P0_mult (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P0.
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
  ! DIMENSION(3,npoints)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
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

  pure subroutine elem_P0_sim (celement, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>

  ! Element type identifier. Must be =EL_P0.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelemens)
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
  ! DIMENSION(3,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
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
! Element subroutines for parametric P1 element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_P1 (celement, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
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

  ! Barycentric coordinates of the point where to evaluate
  real(DP), dimension(3), intent(in) :: Dpoint
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

  real(DP) :: dxj !auxiliary variable

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
  ! That is even faster than when using three IF commands for preventing
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

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P1_mult (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

  !<input>

  ! Element type identifier. Must be =EL_P1.
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
  ! DIMENSION(3,npoints).
  !  Dpoints(1,.)=1st barycentric coordinate
  !  Dpoints(2,.)=2nd barycentric coordinate
  !  Dpoints(3,.)=3rd barycentric coordinate
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

  real(DP),dimension(npoints) :: dxj ! auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC,i) = Dpoints(1,i)
      Dbas(2,DER_FUNC,i) = Dpoints(2,i)
      Dbas(3,DER_FUNC,i) = Dpoints(3,i)
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    dxj = 1E0_DP / Ddetj

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_P1_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, &
                          Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1.
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
  ! DIMENSION(3,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
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

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        Dbas(1,DER_FUNC,i,j) = Dpoints(1,i,j)
        Dbas(2,DER_FUNC,i,j) = Dpoints(2,i,j)
        Dbas(3,DER_FUNC,i,j) = Dpoints(3,i,j)
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1E0_DP / Ddetj(1:npoints,j)

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = -(Djac(4,i,j)-Djac(2,i,j))*dxj(i)
          Dbas(2,DER_DERIV_X,i,j) =  Djac(4,i,j)*dxj(i)
          Dbas(3,DER_DERIV_X,i,j) = -Djac(2,i,j)*dxj(i)
!        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
!        do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = (Djac(3,i,j)-Djac(1,i,j))*dxj(i)
          Dbas(2,DER_DERIV_Y,i,j) = -Djac(3,i,j)*dxj(i)
          Dbas(3,DER_DERIV_Y,i,j) =  Djac(1,i,j)*dxj(i)
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric P2 element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_P2 (celement, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
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

  ! Barycentric coordinates of the point where to evaluate
  real(DP), dimension(3), intent(in) :: Dpoint
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

  real(DP) :: dxj,dp1,dp2,dp3 !auxiliary variables

  ! The P2_2D element is specified by six polynomials on the reference element.
  ! These six polynomials are:
  !
  ! p1(xi1,xi2,xi3) = xi1 * (2 * xi1 - 1)
  ! p2(xi1,xi2,xi3) = xi2 * (2 * xi2 - 1)
  ! p3(xi1,xi2,xi3) = xi3 * (2 * xi3 - 1)
  ! p4(xi1,xi2,xi3) = 4 * xi1 * xi2
  ! p5(xi1,xi2,xi3) = 4 * xi2 * xi3
  ! p6(xi1,xi2,xi3) = 4 * xi1 * xi3

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P2-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
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

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P2_mult (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

  !<input>

  ! Element type identifier. Must be =EL_P2.
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
  ! DIMENSION(3,npoints).
  !  Dpoints(1,.)=1st barycentric coordinate
  !  Dpoints(2,.)=2nd barycentric coordinate
  !  Dpoints(3,.)=3rd barycentric coordinate
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

  real(DP),dimension(npoints) :: Dxj,Dp1,Dp2,Dp3 ! auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P2-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  if (Bder(DER_FUNC)) then
    do i=1,npoints
      Dbas(1,DER_FUNC,i)= Dpoints(1,i)*(Dpoints(1,i)-Dpoints(2,i)-Dpoints(3,i))
      Dbas(2,DER_FUNC,i)= Dpoints(2,i)*(Dpoints(2,i)-Dpoints(1,i)-Dpoints(3,i))
      Dbas(3,DER_FUNC,i)= Dpoints(3,i)*(Dpoints(3,i)-Dpoints(1,i)-Dpoints(2,i))
      Dbas(4,DER_FUNC,i)= 4.0_DP*Dpoints(1,i)*Dpoints(2,i)
      Dbas(5,DER_FUNC,i)= 4.0_DP*Dpoints(2,i)*Dpoints(3,i)
      Dbas(6,DER_FUNC,i)= 4.0_DP*Dpoints(1,i)*Dpoints(3,i)
    end do
  endif

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    Dxj(:) = 1E0_DP / Ddetj(1:npoints)
    do i=1,npoints
      Dp1(i)=(1.0_DP-4.0_DP*Dpoints(1,i))*Dxj(i)
      Dp2(i)=(4.0_DP*Dpoints(2,i)-1.0_DP)*Dxj(i)
      Dp3(i)=(4.0_DP*Dpoints(3,i)-1.0_DP)*Dxj(i)
    end do

    !x-derivatives on current element
    if (Bder(DER_DERIV_X)) then
      do i=1,npoints
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
      end do
    endif

    !y-derivatives on current element
    if (Bder(DER_DERIV_Y)) then
      do i=1,npoints
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
      end do
    endif
  endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_P2_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, &
                          Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P2.
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
  ! DIMENSION(3,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
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

  real(DP),dimension(npoints) :: dxj,Dp1,Dp2,Dp3 !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
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
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj,dp1,dp2,dp3) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1E0_DP / Ddetj(1:npoints,j)

      do i=1,npoints
        Dp1(i)=(1.0_DP-4.0_DP*Dpoints(1,i,j))*Dxj(i)
        Dp2(i)=(4.0_DP*Dpoints(2,i,j)-1.0_DP)*Dxj(i)
        Dp3(i)=(4.0_DP*Dpoints(3,i,j)-1.0_DP)*Dxj(i)
      end do


      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
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
!        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
!        do i=1,npoints
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
        end do
!      ENDIF

    end do
    !$omp end parallel do

  end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric P1~ element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_P1T (celement, Dcoords, Djac, ddetj, Bder, &
                            Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1T.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
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

  ! Barycentric coordinates of the point where to evaluate
  real(DP), dimension(3), intent(in) :: Dpoint
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

  real(DP) :: dxj !auxiliary variable

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
  ! That is even faster than when using three IF commands for preventing
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
      Dbas(3,DER_DERIV_X) = -2E0_DP*Djac(4)*dxj
!    endif

    !y-derivatives on current element
!    if (el_bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -2E0_DP*Djac(1)*dxj
      Dbas(2,DER_DERIV_Y) = -2E0_DP*(Djac(3)- Djac(1))*dxj
      Dbas(3,DER_DERIV_Y) =  2E0_DP*Djac(3)*dxj
!    endif
!  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P1T_mult (celement, Dcoords, Djac, Ddetj, &
                                 Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

  !<input>

  ! Element type identifier. Must be =EL_P1T.
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
  ! DIMENSION(3,npoints).
  !  Dpoints(1,.)=1st barycentric coordinate
  !  Dpoints(2,.)=2nd barycentric coordinate
  !  Dpoints(3,.)=3rd barycentric coordinate
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

  real(DP),dimension(npoints) :: dxj ! auxiliary variable

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1~-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC,i) = 1._DP -2E0_DP*Dpoints(3,i)
      Dbas(2,DER_FUNC,i) = 1._DP -2E0_DP*Dpoints(1,i)
      Dbas(3,DER_FUNC,i) = 1._DP -2E0_DP*Dpoints(2,i)
    end do
  !ENDIF

  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj(:) = 1E0_DP / Ddetj(1:npoints)

    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
        Dbas(1,DER_DERIV_X,i) =  2E0_DP*Djac(2,i)*dxj(i)
        Dbas(2,DER_DERIV_X,i) =  2E0_DP*(Djac(4,i)-Djac(2,i))*dxj(i)
        Dbas(3,DER_DERIV_X,i) = -2E0_DP*Djac(4,i)*dxj(i)
!      END DO
!    ENDIF

    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = -2E0_DP*Djac(1,i)*dxj(i)
        Dbas(2,DER_DERIV_Y,i) = -2E0_DP*(Djac(3,i)- Djac(1,i))*dxj(i)
        Dbas(3,DER_DERIV_Y,i) =  2E0_DP*Djac(3,i)*dxj(i)
      end do
!    ENDIF
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_P1T_sim (celement, Dcoords, Djac, Ddetj, &
                           Bder, Dbas, npoints, nelements, &
                           Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1T.
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
  ! DIMENSION(3,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
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

  real(DP),dimension(npoints) :: dxj !auxiliary variable

  integer :: i   ! point counter
  integer :: j   ! element counter

  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  if (Bder(DER_FUNC)) then

    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements

      do i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 1._DP -2E0_DP*Dpoints(3,i,j)
        Dbas(2,DER_FUNC,i,j) = 1._DP -2E0_DP*Dpoints(1,i,j)
        Dbas(3,DER_FUNC,i,j) = 1._DP -2E0_DP*Dpoints(2,i,j)
      end do

    end do
    !$omp end parallel do

  end if

  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

    !$omp parallel do default(shared) private(i,dxj) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      Dxj(:) = 1E0_DP / Ddetj(1:npoints,j)

      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        do i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) =  2E0_DP*Djac(2,i,j)*Dxj(i)
          Dbas(2,DER_DERIV_X,i,j) =  2E0_DP*(Djac(4,i,j)-Djac(2,i,j))*Dxj(i)
          Dbas(3,DER_DERIV_X,i,j) = -2E0_DP*Djac(4,i,j)*Dxj(i)
!        end do
!      ENDIF

      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
!        do i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -2E0_DP*Djac(1,i,j)*Dxj(i)
          Dbas(2,DER_DERIV_Y,i,j) = -2E0_DP*(Djac(3,i,j)- Djac(1,i,j))*Dxj(i)
          Dbas(3,DER_DERIV_Y,i,j) =  2E0_DP*Djac(3,i,j)*Dxj(i)
        end do
!      ENDIF

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
  ! Lowest order Raviart-Thomas element
  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_RT1_2D (celement, reval, Bder, Dbas)

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
  !
  ! Due to the fact that this is a vector valued basis function, the
  ! meaning of Dbas is extended. There is
  !  Dbas(i        ,:,:,:) = values of the first basis function
  !  Dbas(i+ndofloc,:,:,:) = values of the 2nd basis function, 
  ! with ndofloc the number of local DOFs per element.
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The lowest order Raviart-Thomas element is a vector-valued finite element.
  ! On each cell, one has two P1 polynoms defined as follows:
  !
  !   v = ( v_1 ) = ( a + c x )
  !       ( v_2 )   ( b + c y )
  !
  ! with a,b,c three coefficients. The element has three degrees of
  ! freedom per element which are, however, not a,b,c but the normal
  ! flux through each edge of the triangle:
  !
  !      x3
  !      |\
  !      | \  _
  !      |  \ /|
  !    <-3   2
  !      |    \
  !      |     \
  !      |      \
  !     x1---1---x2
  !          |
  !          v
  !
  ! The three DOFs are 
  !    d1=<v,n1>, 
  !    d2=<v,n2> and 
  !    d3=<v,n3> with
  ! ni the three normal vectors.
  !
  ! We apply a nonparametric approach here...
  
    ! Local variables
    real(DP), dimension(3,3) :: Da,Dcoeff
    real(DP) :: dlen1,dlen2,dlen3
    real(DP) :: dx1,dy1,dx2,dy2,dx3,dy3, dx,dy
    real(DP) :: dnx1,dny1,dnx2,dny2,dnx3,dny3
    real(DP) :: dmx1,dmy1,dmx2,dmy2,dmx3,dmy3
    real(DP) :: da1,da2,da3,db1,db2,db3,dc1,dc2,dc3
    integer :: itwist,i,iel
    real(DP) :: dtwist1,dtwist2,dtwist3
    logical :: bsuccess

    ! Loop through all elements
    !$omp parallel do default(shared)&
    !$omp private(dx1,dx2,dx3,dy1,dy2,dy3,dlen1,dlen2,dlen3,&
    !$omp         dnx1,dnx2,dnx3,dny1,dny2,dny3,dmx1,dmx2,dmx3,&
    !$omp         dmy1,dmy2,dmy3,itwist,dtwist1,dtwist2,dtwist3)&
    !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
    do iel = 1, reval%nelements

      ! Get the twist indices that define the orientation of our edges.
      ! A value of 1 is standard, a value of -1 results in a change of the
      ! sign in the basis functions.
      ! itwistIndex is a bitfield. Each bit specifies the orientation of an edge.
      ! We use the bit to calculate a "1.0" if the edge has positive orientation
      ! and "-1.0" if it has negative orientation.
      itwist = reval%p_ItwistIndex(iel)
      dtwist1 = real(1-iand(int(ishft(itwist, 1)),2),DP)
      dtwist2 = real(1-iand(int(ishft(itwist, 0)),2),DP)
      dtwist3 = real(1-iand(int(ishft(itwist,-1)),2),DP)

      ! Calculate the three normal vectors...
      !
      ! Element corners
      dx1 = reval%p_Dcoords(1,1,iel)
      dy1 = reval%p_Dcoords(2,1,iel)
      dx2 = reval%p_Dcoords(1,2,iel)
      dy2 = reval%p_Dcoords(2,2,iel)
      dx3 = reval%p_Dcoords(1,3,iel)
      dy3 = reval%p_Dcoords(2,3,iel)
      
      ! Inverse of the edge length
      dlen1 = 1.0_DP/sqrt( (dx2-dx1)**2 + (dy2-dy1)**2 )
      dlen2 = 1.0_DP/sqrt( (dx3-dx2)**2 + (dy3-dy2)**2 )
      dlen3 = 1.0_DP/sqrt( (dx1-dx3)**2 + (dy1-dy3)**2 )
      
      ! Calculate the normal vectors to the edges from the corners.
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
      dnx1 = (dy2-dy1)  *dlen1*dtwist1
      dny1 = -(dx2-dx1) *dlen1*dtwist1

      dnx2 = (dy3-dy2)  *dlen2*dtwist2
      dny2 = -(dx3-dx2) *dlen2*dtwist2

      dnx3 = (dy1-dy3)  *dlen3*dtwist3
      dny3 = -(dx1-dx3) *dlen3*dtwist3
      
      ! Midpoints of the edges
      dmx1 = 0.5_DP*(dx1+dx2)
      dmx2 = 0.5_DP*(dx2+dx3)
      dmx3 = 0.5_DP*(dx3+dx1)

      dmy1 = 0.5_DP*(dy1+dy2)
      dmy2 = 0.5_DP*(dy2+dy3)
      dmy3 = 0.5_DP*(dy3+dy1)
      
      ! Set up a transformation matrix to calculate the a,b and c.
      !
      ! The conditions to be solved for one basis function read:
      !
      !   < v,ni > = di
      !
      ! or more exactly
      !
      !   n1x (a+cx)  +  n1y (b+cy)  =  d1
      !   n2x (a+cx)  +  n2y (b+cy)  =  d2
      !   n3x (a+cx)  +  n3y (b+cy)  =  d3
      !
      ! or equivalently
      !
      !  ( n1x  n1y  (n1x x + n1y y) )  ( a )  =  ( d1 )
      !  ( n2x  n2y  (n2x x + n2y y) )  ( b )     ( d2 )
      !  ( n3x  n3y  (n3x x + n3y y) )  ( c )     ( d3 )
      !
      ! Set up this matrix. For (x,y), insert the midpoint of the edge
      ! corresponding to the normal.
      
      Da(1,1) = dnx1
      Da(2,1) = dnx2
      Da(3,1) = dnx3

      Da(1,2) = dny1
      Da(2,2) = dny2
      Da(3,2) = dny3

      Da(1,3) = dnx1*dmx1 + dny1*dmy1
      Da(2,3) = dnx2*dmx2 + dny2*dmy2
      Da(3,3) = dnx3*dmx3 + dny3*dmy3
      
      ! Taking the inverse will give us all a,b and c of
      ! the three basis functions:
      !
      !  ( n1x  n1y  (n1x x + n1y y) )  ( a1  a2  a3 )  =  ( 1  0  0 )
      !  ( n2x  n2y  (n2x x + n2y y) )  ( b1  b2  b3 )     ( 0  1  0 )
      !  ( n3x  n3y  (n3x x + n3y y) )  ( c1  c2  c3 )     ( 0  0  1 )
      !
      !                                 --------------
      !                                 Inverse of the matrix
      call mprim_invert3x3MatrixDirect(Da,Dcoeff,bsuccess)
      
      da1 = Dcoeff(1,1)
      db1 = Dcoeff(2,1)
      dc1 = Dcoeff(3,1)

      da2 = Dcoeff(1,2)
      db2 = Dcoeff(2,2)
      dc2 = Dcoeff(3,2)

      da3 = Dcoeff(1,3)
      db3 = Dcoeff(2,3)
      dc3 = Dcoeff(3,3)

      ! Loop through all points on the current element
      do i = 1, reval%npointsPerElement

        ! Get the point coordinates
        dx = reval%p_DpointsReal(1,i,iel)
        dy = reval%p_DpointsReal(2,i,iel)

        ! Evaluate basis functions.
        !
        !   v = ( a + c*x )
        !       ( b + c*y )
        Dbas(1,DER_FUNC2D,i,iel) = da1 + dc1*dx
        Dbas(2,DER_FUNC2D,i,iel) = da2 + dc2*dx
        Dbas(3,DER_FUNC2D,i,iel) = da3 + dc3*dx

        Dbas(4,DER_FUNC2D,i,iel) = db1 + dc1*dy
        Dbas(5,DER_FUNC2D,i,iel) = db2 + dc2*dy
        Dbas(6,DER_FUNC2D,i,iel) = db3 + dc3*dy
        
        ! 1st derivative is also directly available and
        ! costs no additional time to be computed.
        Dbas(1,DER_DERIV2D_X,i,iel) = dc1
        Dbas(2,DER_DERIV2D_X,i,iel) = dc2
        Dbas(3,DER_DERIV2D_X,i,iel) = dc3

        Dbas(4,DER_DERIV2D_X,i,iel) = 0.0_DP
        Dbas(5,DER_DERIV2D_X,i,iel) = 0.0_DP
        Dbas(6,DER_DERIV2D_X,i,iel) = 0.0_DP

        Dbas(1,DER_DERIV2D_Y,i,iel) = 0.0_DP
        Dbas(2,DER_DERIV2D_Y,i,iel) = 0.0_DP
        Dbas(3,DER_DERIV2D_Y,i,iel) = 0.0_DP

        Dbas(4,DER_DERIV2D_Y,i,iel) = dc1
        Dbas(5,DER_DERIV2D_Y,i,iel) = dc2
        Dbas(6,DER_DERIV2D_Y,i,iel) = dc3

      end do ! i

    end do ! iel
    !$omp end parallel do
     
  ! From these three DOFs, one can extract the three coefficients
  ! a,b,c that define the polynoms of the x- and y-velocity as
  ! follows. One takes into account that the d_i are constant
  ! along the complete edge, including the corner.
  ! So there is for the reference element
  !
  !    d1 = <v,n1> = -(b+cy)
  !    d2 = <v,n2> = 1/sqrt(2) ( a + b + cx + cy )
  !    d3 = <v,n3> = -(a+cx)
  ! 
  ! On the reference triangle, this implies
  !
  !    d1 = -b                        due to y=0
  !    d2 = 1/sqrt(2) ( a + b + c )   for x=y=1/2 on the edge
  !    d3 = -a                        due to x=0
  !
  ! and thus,
  !
  !    b = -d1
  !    a = -d3
  !    c = sqrt(2)*d2 + d1 + d3
  !
  ! The coefficients of the three local basis functions per element 
  ! therefore read
  !
  !   d1=1:  a1 = 0 ,  b1 = -1,  c1 = 1
  !   d2=1:  a2 = 0 ,  b2 = 0 ,  c2 = sqrt(2)
  !   d3=1:  s3 = -1,  b3 = 0 ,  c3 = 1
  !
  ! However, one problem remains - the coordinates are given in
  ! barycentric coordinates. However, this is not really a
  ! problem since
  !
  !    xi1 = 1-x-y
  !    xi2 = x
  !    xi3 = y
  !
  ! so we can just work with xi1 and xi2 instead of x and y.
  !
  ! WARNING: The DOF corresponds to the amount of flow in normal
  ! direction of the edge. However, "normal direction" refers
  ! to the *global* orientation of the edge. Hence, we have to take
  ! this global orientation into account, which is realised by
  ! the "twist indices". If the bit in the twist index corresponding
  ! to the edge is =0, the basis function to this edge is defined as
  ! above. if it is =-1, the orientation is negative and we have
  ! to switch the sign in the basis functions appropriately!
  !
  ! The construction that takes the twist indices into account
  ! reads as follows. For the reference element, there is
  !
  !    d1 = <v,d1 n1> = -(b+cy) * t1
  !    d2 = <v,d2 n2> = 1/sqrt(2) ( a + b + cx + cy ) * t2
  !    d3 = <v,d3 n3> = -(a+cx) * t3
  !
  ! with ti=1 or =-1 depending on the orientation of edge i.
  ! 
  ! On the reference triangle, this implies
  !
  !    d1 = -b * t1                       due to y=0
  !    d2 = 1/sqrt(2) ( a + b + c ) * t2  for x=y=1/2 on the edge
  !    d3 = -a * t3                       due to x=0
  !
  ! and thus,
  !
  !    b = -d1 * t1
  !    a = -d3 * t3
  !    c = sqrt(2)*d2*t2 + d1*t1 + d3*t3
  !
  ! The coefficients of the three local basis functions per element 
  ! therefore read
  !
  !   d1=1:  a1 = 0  ,  b1 = -t1,  c1 = t1
  !   d2=1:  a2 = 0  ,  b2 = 0  ,  c2 = sqrt(2)*t2
  !   d3=1:  a3 = -t3,  b3 = 0  ,  c3 = t3
  !

!  ! Local variables
!  real(DP) :: dxi1,dxi2,dxj,dtwist1,dtwist2,dtwist3
!  integer :: i,j
!  integer :: itwist
!
!    ! Calculate function values?
!    if(Bder(DER_FUNC2D)) then
!
!      ! If function values are desired, calculate them.
!
!      ! Loop through all elements
!      !$omp parallel do default(shared) private(i,dxi1,dxi2,itwist,dtwist1,dtwist2,dtwist3)&
!      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
!      do j = 1, reval%nelements
!
!        ! Get the twist indices that define the orientation of our edges.
!        ! A value of 1 is standard, a value of -1 results in a change of the
!        ! sign in the basis functions.
!        ! itwistIndex is a bitfield. Each bit specifies the orientation of an edge.
!        ! We use the bit to calculate a "1.0" if the edge has positive orientation
!        ! and "-1.0" if it has negative orientation.
!        itwist = reval%p_ItwistIndex(j)
!        dtwist1 = real(1-iand(int(ishft(itwist, 1)),2),DP)
!        dtwist2 = real(1-iand(int(ishft(itwist, 0)),2),DP)
!        dtwist3 = real(1-iand(int(ishft(itwist,-1)),2),DP)
!
!        ! Loop through all points on the current element
!        do i = 1, reval%npointsPerElement
!
!          ! Get the point coordinates
!          dxi1 = reval%p_DpointsRef(2,i,j)
!          dxi2 = reval%p_DpointsRef(3,i,j)
!
!          ! Evaluate basis functions
!          !
!          !   v1 = (       t1 x )
!          !        ( -t1 + t1 y )
!          !
!          !   v2 = ( sqrt(2) t2 x )
!          !        ( sqrt(2) t2 y )
!          !
!          !   v3 = ( -t3 + t3 x )
!          !        (       t3 y )
!          Dbas(1,DER_FUNC2D,i,j) = dtwist1*dxi1
!          Dbas(2,DER_FUNC2D,i,j) = dtwist2*sqrt(2.0_DP)*dxi1
!          Dbas(3,DER_FUNC2D,i,j) = dtwist3*(-1.0_DP + dxi1)
!
!          Dbas(4,DER_FUNC2D,i,j) = dtwist1*(-1.0_DP + dxi2)
!          Dbas(5,DER_FUNC2D,i,j) = dtwist2*sqrt(2.0_DP)*dxi2
!          Dbas(6,DER_FUNC2D,i,j) = dtwist3*dxi2
!
!        end do ! i
!
!      end do ! j
!      !$omp end parallel do
!     
!    end if
!
!    ! Calculate derivatives?
!    if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
!
!      ! If function values are desired, calculate them.
!
!      ! Loop through all elements
!      !$omp parallel do default(shared) private(i,dxi1,dxi2,dxj)&
!      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
!      do j = 1, reval%nelements
!
!        ! Loop through all points on the current element
!        do i = 1, reval%npointsPerElement
!
!          ! Get the point coordinates
!          dxi1 = reval%p_DpointsRef(1,i,j)
!          dxi2 = reval%p_DpointsRef(2,i,j)
!          
!          ! Multiplier for the inverse of the mapping
!          dxj = 1.0_DP / reval%p_Ddetj(i,j)
!
!          ! Evaluate basis functions
!          !
!          ! v1 = (0,-1) + (x,y)
!          ! v2 = sqrt(2)* (x,y)
!          ! v3 = (-1,0) + (x,y)
!          !
!          ! => Dv1 = Dv3 = (1 0) , Dv2 = ( sqrt(2)          )
!          !                (0 1)         (          sqrt(2) )
!          !
!          ! Dvi is multiplied to the mapping matrix which gives the
!          ! actual derivative:
!          !    (1 0) * (e f) = (e f) = 1/(ad-bc) (  d -b )
!          !    (0 1)   (g h)   (g h)             ( -c  a )
!          Dbas(1,DER_DERIV_X,i,j) =  dtwist1*reval%p_Djac(4,i,j)*dxj
!          Dbas(2,DER_DERIV_X,i,j) =  dtwist2*sqrt(2.0_DP)*reval%p_Djac(4,i,j)*dxj
!          Dbas(3,DER_DERIV_X,i,j) =  dtwist3*reval%p_Djac(4,i,j)*dxj
!
!          Dbas(4,DER_DERIV_X,i,j) = -dtwist1*reval%p_Djac(2,i,j)*dxj
!          Dbas(5,DER_DERIV_X,i,j) = -dtwist2*sqrt(2.0_DP)*reval%p_Djac(2,i,j)*dxj
!          Dbas(6,DER_DERIV_X,i,j) = -dtwist3*reval%p_Djac(2,i,j)*dxj
!
!          Dbas(1,DER_DERIV_Y,i,j) = -dtwist1*reval%p_Djac(3,i,j)*dxj
!          Dbas(2,DER_DERIV_Y,i,j) = -dtwist2*sqrt(2.0_DP)*reval%p_Djac(3,i,j)*dxj
!          Dbas(3,DER_DERIV_Y,i,j) = -dtwist3*reval%p_Djac(3,i,j)*dxj
!
!          Dbas(4,DER_DERIV_Y,i,j) =  dtwist1*reval%p_Djac(1,i,j)*dxj
!          Dbas(5,DER_DERIV_Y,i,j) =  dtwist2*sqrt(2.0_DP)*reval%p_Djac(1,i,j)*dxj
!          Dbas(6,DER_DERIV_Y,i,j) =  dtwist3*reval%p_Djac(1,i,j)*dxj
!
!        end do ! i
!
!      end do ! j
!      !$omp end parallel do


  end subroutine

  !************************************************************************
  ! Discontinuous P1 element
  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

 subroutine elem_eval_DCP1_2D (celement, reval, Bder, Dbas)

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
  !
  ! Due to the fact that this is a vector valued basis function, the
  ! meaning of Dbas is extended. There is
  !  Dbas(i        ,:,:,:) = values of the first basis function
  !  Dbas(i+ndofloc,:,:,:) = values of the 2nd basis function, 
  ! with ndofloc the number of local DOFs per element.
  real(DP), dimension(:,:,:,:), intent(out)      :: Dbas
!</output>

! </subroutine>

    real(DP), dimension(reval%npointsPerElement) :: dxj !auxiliary variable

    integer :: i   ! point counter
    integer :: j   ! element counter

    ! Clear the output array
    !Dbas = 0.0_DP

    !if function values are desired
    if (Bder(DER_FUNC)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,reval%nelements
      
        ! There are three basis functions on the reference element:
        !   phi_1 = 1
        !   phi_2 = xi   = lambda_2 (in barycentric coordinates)
        !   phi_3 = eta  = lambda_3 (in barycentric coordinates)

        do i=1,reval%npointsPerElement
          Dbas(1,DER_FUNC,i,j) = 1.0_DP
          Dbas(2,DER_FUNC,i,j) = reval%p_DpointsRef(2,i,j)
          Dbas(3,DER_FUNC,i,j) = reval%p_DpointsRef(3,i,j)
        end do

      end do
      !$omp end parallel do

    end if

    !if x-or y-derivatives are desired
    if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then

      !$omp parallel do default(shared) private(i,dxj) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,reval%nelements
        Dxj(:) = 1E0_DP / reval%p_Ddetj(1:reval%npointsPerElement,j)

        !x-derivatives on current element
  !      IF (Bder(DER_DERIV_X)) THEN
          do i=1,reval%npointsPerElement
            Dbas(1,DER_DERIV_X,i,j) = 0.0_DP
            Dbas(2,DER_DERIV_X,i,j) =  reval%p_Djac(4,i,j)*dxj(i)
            Dbas(3,DER_DERIV_X,i,j) = -reval%p_Djac(2,i,j)*dxj(i)
  !        end do
  !      ENDIF

        !y-derivatives on current element
  !      IF (Bder(DER_DERIV_Y)) THEN
  !        do i=1,npoints
            Dbas(1,DER_DERIV_Y,i,j) = 0.0_DP
            Dbas(2,DER_DERIV_Y,i,j) = -reval%p_Djac(3,i,j)*dxj(i)
            Dbas(3,DER_DERIV_Y,i,j) =  reval%p_Djac(1,i,j)*dxj(i)
          end do
  !      ENDIF

      end do
      !$omp end parallel do

    end if
     

  end subroutine

end module
