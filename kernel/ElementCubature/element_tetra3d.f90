!##############################################################################
!# ****************************************************************************
!# <name> element_tetra3d </name>
!# ****************************************************************************
!# 
!# <purpose>
!# This module contains the implementations of the 3D tetrahedron basis
!# functions.
!#
!# </purpose>
!##############################################################################

module element_tetra3d

  use fsystem
  use elementbase
  use derivatives

  implicit none
  
  private
  
  public :: elem_P0_3D 
  public :: elem_P0_3D_mult 
  public :: elem_P0_3D_sim 
  public :: elem_P1_3D 
  public :: elem_P1_3D_mult 
  public :: elem_P1_3D_sim 

contains

!**************************************************************************
! Element subroutines for parametric P0_3D element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  pure subroutine elem_P0_3D (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P0_3D.
  integer(I32), intent(in)  :: celement
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
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
  real(DP), dimension(4), intent(in) :: Dpoint
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

  pure subroutine elem_P0_3D_mult (celement, Dcoords, Djac, Ddetj, &
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
  ! DIMENSION(4,npoints)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
  !  Dpoints(4,.) = 4th barycentric coordinate
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

  pure subroutine elem_P0_3D_sim (celement, Dcoords, Djac, Ddetj, &
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
  real(DP), dimension(:,:), intent(in) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
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

!</subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  end subroutine 

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
! Let us use the notation 
!
!    P(X,Y,Z)  :=  p( t^-1 (X,Y,Z) )  =  p (xi1,xi2,xi3,xi4)
!
! which is the 'standard' form of the polynomials on the
! reference element without barycentric coordinates involved (i.e.
! for the P1-element it is of the form P(x,y,z) = c1*x + c2*y + c3*z + c4).
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

  pure subroutine elem_P1_3D (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1_3D.
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
  !  Djac(3,i) = J_i(3,1)
  !  Djac(4,i) = J_i(1,2)
  !  Djac(5,i) = J_i(2,2)
  !  Djac(6,i) = J_i(3,2)
  !  Djac(7,i) = J_i(1,3)
  !  Djac(8,i) = J_i(2,3)
  !  Djac(9,i) = J_i(3,3)
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
  real(DP), dimension(4), intent(in) :: Dpoint
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
  
  ! A matrix for the inverse jacobian matrix
  real(DP), dimension(9) :: Dinv
  
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
  ! That is even faster than when using three IF commands for preventing
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
    
  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_P1_3D_mult (celement, Dcoords, Djac, Ddetj, &
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
  !  Djac(3,i) = J_i(3,1)
  !  Djac(4,i) = J_i(1,2)
  !  Djac(5,i) = J_i(2,2)
  !  Djac(6,i) = J_i(3,2)
  !  Djac(7,i) = J_i(1,3)
  !  Djac(8,i) = J_i(2,3)
  !  Djac(9,i) = J_i(3,3)
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

  real(DP) :: dxj ! auxiliary variable
  
  integer :: i   ! point counter

  ! A matrix for the inverse jacobian matrix
  real(DP), dimension(9) :: Dinv
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC3D,i) = Dpoints(1,i)
      Dbas(2,DER_FUNC3D,i) = Dpoints(2,i)
      Dbas(3,DER_FUNC3D,i) = Dpoints(3,i)
      Dbas(4,DER_FUNC3D,i) = Dpoints(4,i)
    end do
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
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF
    
  end subroutine 

  !************************************************************************
  
!<subroutine>  

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_P1_3D_sim (celement, Dcoords, Djac, Ddetj, &
                             Bder, Dbas, npoints, nelements, Dpoints)

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

  real(DP) :: dxj !auxiliary variable
  
  integer :: i   ! point counter
  integer :: j   ! element counter

  ! A matrix for the inverse jacobian matrix
  real(DP), dimension(9) :: Dinv
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  if (Bder(DER_FUNC3D)) then
    
    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > EL_NELEMMIN_OMP)
    do j=1,nelements
    
      do i=1,npoints
        Dbas(1,DER_FUNC3D,i,j) = Dpoints(1,i,j)
        Dbas(2,DER_FUNC3D,i,j) = Dpoints(2,i,j)
        Dbas(3,DER_FUNC3D,i,j) = Dpoints(3,i,j)
        Dbas(4,DER_FUNC3D,i,j) = Dpoints(4,i,j)
      end do
      
    end do
    !$omp end parallel do
    
  end if
    
  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
      (Bder(DER_DERIV3D_Z))) then
  
    !$omp parallel do default(shared) private(i,dxj,Dinv) &
    !$omp if(nelements > EL_NELEMMIN_OMP)
    do j=1,nelements

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
        do i=1,npoints
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
        end do
!      ENDIF
    end do
    !$omp end parallel do
      
  end if
    
  end subroutine 

end module
