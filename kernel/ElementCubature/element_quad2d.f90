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

  use elementbase
  use derivatives
  use mprimitives
  use transformation

implicit none

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

  pure subroutine elem_Q0 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q0.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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

  pure subroutine elem_Q0_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
  
  !</input>
  
  !<output>
  
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
  
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

  pure subroutine elem_Q0_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
  
  !</input>
  
  !<output>
  
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
  
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

  pure subroutine elem_Q1 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
    
  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_Q1_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable
  
  integer :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
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
    dxj = 0.25E0_DP / Ddetj
    
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

  pure subroutine elem_Q1_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC)) then
  
    !$omp parallel do default(shared) private(i)
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
  
    !$omp parallel do default(shared) private(i,dxj,Dhelp)
    do j=1,nelements
      dxj = 0.25E0_DP / Ddetj(:,j)
      
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

  pure subroutine elem_Q2 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
    ! Let's assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It's inverse s^{-1} maps to the reference element.
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
    ! This is computable now. Let's write:
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

  pure subroutine elem_Q2_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
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
    Dxj = 1.0E0_DP / Ddetj
    
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
    ! Let's assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It's inverse s^{-1} maps to the reference element.
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
    ! This is computable now. Let's write:
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
    Dxj = 1.0E0_DP / Ddetj
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

  pure subroutine elem_Q2_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC)) then
  
    !$omp parallel do default(shared) private(i,dx,dy)
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
  
    !$omp parallel do default(shared) private(i,Dxj,dx,dy,Dhelp,idof)
    do j=1,nelements
      Dxj = 1.0E0_DP / Ddetj(:,j)
      
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
    ! Let's assume, out finite element function on the real element is
    ! f(x) and the corresponding function on the reference element g(y),
    ! so x is the real coordinate and y the reference coordinate.
    ! There is a mapping s:[-1,1]->R2 that maps to the real element.
    ! It's inverse s^{-1} maps to the reference element.
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
    ! This is computable now. Let's write:
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
    do j=1,nelements
      Dxj = 1.0E0_DP / Ddetj(:,j)
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

  end if
    
  end subroutine 

!**************************************************************************
! Element subroutines for parametric QP1 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  pure subroutine elem_QP1 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                            Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_QP1.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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

  pure subroutine elem_EM30 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE),
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on the real element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
  !CALL INVERT(A,F,CKH,0)
  call mprim_invertMatrixPivotDble(A,4)

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
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function
    
    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function
    
    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function
    
    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_EM30_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given points. The coordinates are expected
  ! on the real element!
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables  
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
  
  !CALL INVERT(A,F,CKH,0)
  call mprim_invertMatrixPivotDble(A,4)

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
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function
    
    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function
    
    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function
    
    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine 

  !************************************************************************
  
!<subroutine>  

  subroutine elem_EM30_sim (ieltyp, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the real element.
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables  
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
  if (iand(ieltyp,int(2**17,I32)) .eq. 0) then

    ! Check whether to scale the local coordinate system or not.
  
    if (iand(ieltyp,int(2**18,I32)) .eq. 0) then
  
      ! Use pivoting and scaled local coordinate system for 
      ! increased numerical stability.
    
      ! Loop over the elements
      
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
        !CALL INVERT(A,F,CKH,0)
        call mprim_invertMatrixPivotDble(A,4)

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
      
    else
    
      ! Use pivoting for increased numerical stability.
      ! Don't scaled local coordinate system 
    
      ! Loop over the elements
      
      do j=1,nelements
        
        do IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
        end do

        ! Don't scale the local coordinate system in this approach.
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
        !CALL INVERT(A,F,CKH,0)
        call mprim_invertMatrixPivotDble(A,4)

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

    end if
    
  else
  
    ! Don't use pivoting.
    
    ! Check whether to scae the local coordinate system or not.
  
    if (iand(ieltyp,int(2**18,I32)) .eq. 0) then
  
      ! Loop over the elements
      
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
        call mprim_invert4x4MatrixDirectDble(A,CK)

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
    
    else
    
      ! Loop over the elements
      
      do j=1,nelements
        
        do IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,mod(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,mod(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,mod(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,mod(IVE,4)+1,j)-Dcoords(2,IVE,j))
        end do

        ! Don't scale the local coordinate in this approach.
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
        call mprim_invert4x4MatrixDirectDble(A,CK)

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

    end if
  
  end if
             
  contains

    ! Auxiliary functions
      
    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function
    
    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function
    
    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function
    
    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
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

  pure subroutine elem_E030 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
    
  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_E030_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable
  
  integer :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
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
    dxj = 0.125E0_DP / Ddetj
    
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

  pure subroutine elem_E030_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                 Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC)) then
  
    !$omp parallel do default(shared) private(i)
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
  
    !$omp parallel do default(shared) private(i,dxj,Dhelp)
    do j=1,nelements
      dxj = 0.125E0_DP / Ddetj(:,j)
      
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

  pure subroutine elem_EB30 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1TB.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
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

  pure subroutine elem_EB30_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1TB.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
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
    dxj = 0.125E0_DP / Ddetj
    
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

  pure subroutine elem_EB30_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                 Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1TB.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC)) then
  
    !$omp parallel do default(shared) private(i,dx,dy,dxy)
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
  
    !$omp parallel do default(shared) private(i,dxj,dx,dy,Dhelp)
    do j=1,nelements
      dxj = 0.125E0_DP / Ddetj(:,j)
      
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

  pure subroutine elem_EM31 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on the real element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables  
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
  
  ! Calculate the matrix V (=A) with vij = int_ei Fj (x,y) d(x,y).
  ! Loop over the edges.
  do IA = 1,NVE
    ! Set up the coefficients of the linear system to calculate ai, bi, ci and di.
    ! Use the X- and Y-coordinates of the midpointof every edge to evaluate
    ! the Fi.
    A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
  end do
  
  ! Invert that matrix V to get the matrix of the coefficients of the
  ! four polynomials. The matix A (=V) is replaced by the inverse.
  !CALL INVERT(A,F,CKH,0)
  call mprim_invertMatrixPivotDble(A,4)

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
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function
    
    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function
    
    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function
    
    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_EM31_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given points. The coordinates are expected
  ! on the real element!
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables  
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
  
  !CALL INVERT(A,F,CKH,0)
  call mprim_invertMatrixPivotDble(A,4)

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
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function
    
    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function
    
    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function
    
    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine 

  !************************************************************************
  
!<subroutine>  

  subroutine elem_EM31_sim (ieltyp, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the real element.
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  integer, parameter :: NVE = 4

  ! auxiliary variables  
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
  if (iand(ieltyp,int(2**17,I32)) .eq. 0) then
  
    ! Use pivoting for increased numerical stability.
  
    ! Loop over the elements
    
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
      !CALL INVERT(A,F,CKH,0)
      call mprim_invertMatrixPivotDble(A,4)

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
    
  else
  
    ! Don't use pivoting.
    
    ! Loop over the elements
    
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
      call mprim_invert4x4MatrixDirectDble(A,CK)

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
  
  end if
             
  contains

    ! Auxiliary functions
      
    elemental real(DP) function F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    end function
    
    elemental real(DP) function F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    end function
    
    elemental real(DP) function F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    end function
    
    elemental real(DP) function F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    real(DP), intent(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    end function

  end subroutine 

!**************************************************************************
! Element subroutines for parametric Q1~ element, midpoint value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************
 
!<subroutine>  

  pure subroutine elem_E031 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
    
  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_E031_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(4,NDIM2D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: dxj !auxiliary variable
  
  integer :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
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
    dxj = 0.5E0_DP / Ddetj
    
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

  pure subroutine elem_E031_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                 Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC)) then
  
    !$omp parallel do default(shared) private(i)
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
  
    !$omp parallel do default(shared) private(i,dxj,Dhelp)
    do j=1,nelements
      dxj = 0.5E0_DP / Ddetj(:,j)
      
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

  pure subroutine elem_E050 (ieltyp, Dcoords, itwistIndex, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2T.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges. 
  integer(I32), intent(IN) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
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

  pure subroutine elem_E050_mult (ieltyp, Dcoords, itwistIndex, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2T.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges. 
  integer(I32), intent(IN) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
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

  pure subroutine elem_E050_sim (ieltyp, Dcoords, ItwistIndex, Djac, Ddetj, &
                                 Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2T.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! List of twist indices. For every edge/face on every cell, the twist
  ! index defines the orientation of the edge/face.
  ! Array with DIMENSION(1:NVE/NVA,nelements)
  integer(I32), dimension(:), intent(IN) :: ItwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
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
  
    !$omp parallel do default(shared) private(i,dx,dy,d5,d6,d7,d8)
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
  
    !$omp parallel do default(shared) private(i,dxj,Dhelp,dx,dy,d5,d6,d7,d8)
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

  pure subroutine elem_EB50 (ieltyp, Dcoords, itwistIndex, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2TB.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges. 
  integer(I32), intent(IN) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(2), intent(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(OUT) :: Dbas
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
  ! That's even faster than when using three IF commands for preventing
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

  pure subroutine elem_EB50_mult (ieltyp, Dcoords, itwistIndex, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2TB.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
  ! Twist indices bitfield of the element that defines the orientation of
  ! the edges. 
  integer(I32), intent(IN) :: itwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
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

  pure subroutine elem_EB50_sim (ieltyp, Dcoords, ItwistIndex, Djac, Ddetj, &
                                 Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2TB.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
  ! List of twist indices. For every edge/face on every cell, the twist
  ! index defines the orientation of the edge/face.
  ! Array with DIMENSION(1:NVE/NVA,nelements)
  integer(I32), dimension(:), intent(IN) :: ItwistIndex

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,DER_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
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
  
    !$omp parallel do default(shared) private(i,dx,dy,d5,d6,d7,d8)
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
  
    !$omp parallel do default(shared) private(i,dx,dy,d5,d6,d7,d8,dxj,Dhelp)
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


  !****************************************************************************
  !****************************************************************************
  
  ! -------------- NEW ELEMENT INTERFACE IMPLEMENTATIONS FOLLOW --------------
  
  !****************************************************************************
  !****************************************************************************



  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_Q1_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
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
      
      ! Loop through all elements
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
      
    end if
    
    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
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
      
    end if
  
  end subroutine

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_Q2_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
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
      
    end if
    
    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
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
      
    end if
  
  end subroutine

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_Q2H_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
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
      
    end if
    
    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
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
      
    end if
  
  end subroutine

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_QP1_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
!</output>

! </subroutine>

  ! Element Description
  ! -------------------
  ! The QP1_2D element is specified by three polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  ! { 1, x, y }
  !
  ! As the QP1 element is discontinous, the basis polynomials don't have to
  ! fulfill any special conditions - they're simply defined as:
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
      do j = 1, reval%nelements
      
        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement
        
          ! Evaluate basis functions
          Dbas(1,DER_FUNC2D,i,j) = 1.0_DP
          Dbas(2,DER_FUNC2D,i,j) = reval%p_DpointsRef(1,i,j)
          Dbas(3,DER_FUNC2D,i,j) = reval%p_DpointsRef(2,i,j)
        
        end do ! i
      
      end do ! j
      
    end if
    
    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
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
      
    end if
  
  end subroutine

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_E030_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
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
      
    end if
    
    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
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
      
    end if
  
  end subroutine

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_EB30_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
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
      
    end if
    
    ! Calculate derivatives?
    if(Bder(DER_DERIV2D_X) .or. Bder(DER_DERIV2D_Y)) then

      ! Loop through all elements
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
      
    end if
  
  end subroutine
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_EM30_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The EM30_2D element is specified by four polynomials per element.
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
      call mprim_invert4x4MatrixDirectDble(Da, Dc)
      
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
    
    ! That's it

  end subroutine

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_eval_E050_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
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
  ! Numer. Math., 53 (1988), pp. 701�738.
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
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2)
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
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2,ddet,DrefDer)
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

  pure subroutine elem_eval_EB50_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
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
  ! Numer. Math., 53 (1988), pp. 701�738.
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
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2)
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
      !$omp parallel do private(i,itwist,Dtw,dx,dy,dx2,dy2,ddet,DrefDer)
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

  pure subroutine elem_eval_EM50_2D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the
  ! reference element for multiple given elements.
!</description>

!<input>
  ! The element specifier, must be EL_EM50_2D.
  integer(I32), intent(IN)                       :: celement
  
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(IN)             :: reval
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN)              :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(OUT)      :: Dbas
!</output>

!</subroutine>

  ! Element Description
  ! -------------------
  ! The EM50_2D element is specified by nine polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, y, x*y, x^2, y^2, x^2*y, x*y^2, x^3*y - x*y^3 }
  !
  ! see:
  ! J.-P. Hennart, J. Jaffre, and J. E. Roberts;
  ! "A Constructive Method for Deriving Finite Elements of Nodal Type";
  ! Numer. Math., 53 (1988), pp. 701�738.
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
    ! Remark: We'll use the 1D rule to create the 2D rule here
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
      
      ! ...and also calculate the inverse of the element's area.
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

      call mprim_invertMatrixPivotDble(Da, NBAS)
      
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
    
    ! That's it

  end subroutine

end module