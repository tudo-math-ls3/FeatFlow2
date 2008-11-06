!##############################################################################
!# ****************************************************************************
!# <name> element_hexa3d </name>
!# ****************************************************************************
!# 
!# <purpose>
!# This module contains the implementations of the 3D hexahedron basis
!# functions.
!#
!# </purpose>
!##############################################################################

module element_hexa3d

  use elementbase
  use derivatives
  use mprimitives

implicit none

contains

!**************************************************************************
! Element subroutines for parametric 3D Q0 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  pure subroutine elem_Q0_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q0_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  ! Dpoint(3) = z-coordinate
  real(DP), dimension(3), intent(IN) :: Dpoint
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

  pure subroutine elem_Q0_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates,
  ! Dpoints(3,.)=z-coordinates.
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

  pure subroutine elem_Q0_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
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
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
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

  pure subroutine elem_Q1_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate,
  ! Dpoint(3) = z-coordinate
  real(DP), dimension(3), intent(IN) :: Dpoint
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
  real(DP), dimension(8,NDIM3D) :: Dhelp
  real(DP) :: dx,dy,dz, djx, djy, djz

  real(DP) :: dxj !auxiliary variable
  
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
    
  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_Q1_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
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
  real(DP), dimension(8,NDIM3D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: Dxj !auxiliary variable
  real(DP) :: djx, djy, djz
  integer :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    do i=1,npoints
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
    end do
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR.&
!      (Bder(DER_DERIV3D_Z))) THEN
    Dxj = 0.125E0_DP / Ddetj
    
    !x-, y- and z-derivatives on reference element
    do i=1,npoints
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
    end do
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF
    
  end subroutine 

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_Q1_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1_3D.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
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
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(8,NDIM3D,npoints) :: Dhelp
  real(DP) :: djx, djy, djz
  real(DP),dimension(npoints) :: Dxj !auxiliary variable
  
  integer :: i   ! point counter
  integer :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC3D)) then
  
    !$omp parallel do default(shared) private(i)
    do j=1,nelements
    
      do i=1,npoints
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
      end do
      
    end do
    !$omp end parallel do
    
  end if
    
  !if x-, y- or z-derivatives are desired
  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
      (Bder(DER_DERIV3D_Z))) then
  
    !$omp parallel do default(shared) private(i,Dxj,Dhelp,djx,djy,djz)
    do j=1,nelements
      Dxj = 0.125_DP / Ddetj(:,j)
      
      !x-, y- and z-derivatives on reference element
      do i=1,npoints
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
      end do
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV3D_X)) THEN
        do i=1,npoints
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
!        end do
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
!        do i=1,npoints
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
!        end do
!      ENDIF

      !z-derivatives on current element
!      IF (Bder(DER_DERIV3D_Z)) THEN
!        do i=1,npoints
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
        end do
!      ENDIF
    end do
    !$omp end parallel do
      
  end if
    
  end subroutine 

!**************************************************************************
! Element subroutines for parametric Q1~ element, integral mean value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  pure subroutine elem_E030_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                                Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_E030_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate,
  ! Dpoint(3) = z-coordinate
  real(DP), dimension(3), intent(IN) :: Dpoint
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
  real(DP), dimension(6,NDIM3D) :: Dhelp
  ! auxiliary variables
  real(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  real(DP) :: dxj
  real(DP), parameter :: R16 = 1.0_DP/6.0_DP

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
    
  end subroutine 

  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_E030_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1T_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
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
  real(DP), dimension(6,NDIM3D,npoints) :: Dhelp
  ! auxiliary variables
  real(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  real(DP), dimension(npoints) :: Dxj
  real(DP), parameter :: R16 = 1.0_DP/6.0_DP
  integer :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    do i=1,npoints
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
    end do
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR.&
!      (Bder(DER_DERIV3D_Z))) THEN
    Dxj = 1.0_DP / Ddetj
    
    !x-, y- and z-derivatives on reference element
    do i=1,npoints
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
    end do
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF
    
  end subroutine 

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_E030_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1_3D.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
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
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,NDIM3D,npoints) :: Dhelp
  ! auxiliary variables
  real(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  real(DP),dimension(npoints) :: Dxj
  real(DP), parameter :: R16 = 1.0_DP/6.0_DP
  integer :: i   ! point counter
  integer :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC3D)) then
  
    !$omp parallel do default(shared) private(i,dx,dy,dz,dxy,dyz)
    do j=1,nelements
    
      do i=1,npoints
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
      end do
      
    end do
    !$omp end parallel do
    
  end if
    
  !if x-, y- or z-derivatives are desired
  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
      (Bder(DER_DERIV3D_Z))) then
  
    !$omp parallel do default(shared) private(i,dx,dy,dz,djx,djy,djz,Dhelp,dxj)
    do j=1,nelements
      Dxj = 1.0_DP / Ddetj(:,j)
      
      !x-, y- and z-derivatives on reference element
      do i=1,npoints
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
      end do
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV3D_X)) THEN
        do i=1,npoints
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
!        end do
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
!        do i=1,npoints
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
!        end do
!      ENDIF

      !z-derivatives on current element
!      IF (Bder(DER_DERIV3D_Z)) THEN
!        do i=1,npoints
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
        end do
!      ENDIF
    end do
    !$omp end parallel do
      
  end if
    
  end subroutine 

!**************************************************************************
! Element subroutines for parametric Q1~ element, face-midpoint-value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  pure subroutine elem_E031_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                                Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate,
  ! Dpoint(3) = z-coordinate
  real(DP), dimension(3), intent(IN) :: Dpoint
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
  real(DP), dimension(6,NDIM3D) :: Dhelp
  real(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  real(DP), parameter :: R12 = 0.5_DP
  real(DP), parameter :: R13 = 1.0_DP/3.0_DP
  real(DP), parameter :: R16 = 1.0_DP/6.0_DP

  real(DP) :: dxj !auxiliary variable

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
    
  end subroutine 

  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_E031_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1T_3D.
  integer(I32), intent(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
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
  real(DP), dimension(6,NDIM3D,npoints) :: Dhelp

  real(DP),dimension(npoints) :: Dxj !auxiliary variable
  real(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  real(DP), parameter :: R12 = 0.5_DP
  real(DP), parameter :: R13 = 1.0_DP/3.0_DP
  real(DP), parameter :: R16 = 1.0_DP/6.0_DP
  integer :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    do i=1,npoints
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
    end do
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR.&
!      (Bder(DER_DERIV3D_Z))) THEN
    Dxj = 1.0_DP / Ddetj
    
    !x-, y- and z-derivatives on reference element
    do i=1,npoints
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
    end do
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      do i=1,npoints
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
      end do
!    ENDIF
!  ENDIF
    
  end subroutine 

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_E031_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1_3D.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
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
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,NDIM3D,npoints) :: Dhelp
  real(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  real(DP), parameter :: R12 = 0.5_DP
  real(DP), parameter :: R13 = 1.0_DP/3.0_DP
  real(DP), parameter :: R16 = 1.0_DP/6.0_DP
  real(DP),dimension(npoints) :: Dxj !auxiliary variable
  
  integer :: i   ! point counter
  integer :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC3D)) then
  
    !$omp parallel do default(shared) private(i,dx,dy,dz,dxy,dyz)
    do j=1,nelements
    
      do i=1,npoints
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
      end do
      
    end do
    !$omp end parallel do
    
  end if
    
  !if x-, y- or z-derivatives are desired
  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
      (Bder(DER_DERIV3D_Z))) then
  
    !$omp parallel do default(shared) private(i,dx,dy,dz,djx,djy,djz,Dhelp,Dxj)
    do j=1,nelements
      Dxj = 1.0_DP / Ddetj(:,j)
      
      !x-, y- and z-derivatives on reference element
      do i=1,npoints
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
      end do
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV3D_X)) THEN
        do i=1,npoints
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
!        end do
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
!        do i=1,npoints
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
!        end do
!      ENDIF

      !z-derivatives on current element
!      IF (Bder(DER_DERIV3D_Z)) THEN
!        do i=1,npoints
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
        end do
!      ENDIF
    end do
    !$omp end parallel do
      
  end if
    
  end subroutine 

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

  pure subroutine elem_EM30_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                                Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30_3D.
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
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on the real element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  real(DP), dimension(3), intent(IN) :: Dpoint
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

  real(DP), parameter :: Q1 =   0.25_DP
  real(DP), parameter :: Q2 =   1.0_DP /   3.0_DP
  real(DP), parameter :: Q3 =   0.2_DP
  real(DP), parameter :: Q4 =   0.6_DP
  real(DP), parameter :: Q8 =  -9.0_DP /  16.0_DP
  real(DP), parameter :: Q9 = 100.0_DP / 192.0_DP
  
  real(DP), dimension(3,8) :: P                       ! cubature points
  real(DP) :: PP1,PP2,PP3,PP4,PP5,PS1,PS2,PQ1,PQ2,PQ3 ! cubature weights
  real(DP), dimension(6,6) :: A, B                    ! coefficient matrices
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC1,&           ! auxiliary coefficients
              CC2,CC3,CD1,CD2,CD3,CD4,CD5,&
              CD6,CE1,CE2,CE3,CE4,CE5,CE6
  real(DP), dimension(6,10) :: COB                    ! monomial coefficients
  real(DP), dimension(3,4,6) :: V                     ! face corner vertices
  real(DP), dimension(2:6,8) :: F                     ! transformed monomials
  real(DP) :: dx,dy,dz
  integer :: k,l
  
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
    do k=1,6
    
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
      PP1=sqrt((V(1,1,k)-V(1,2,k))**2+(V(2,1,k)-V(2,2,k))**2+&
          (V(3,1,k)-V(3,2,k))**2)
      PP2=sqrt((V(1,2,k)-V(1,3,k))**2+(V(2,2,k)-V(2,3,k))**2+&
          (V(3,2,k)-V(3,3,k))**2)
      PP3=sqrt((V(1,1,k)-V(1,3,k))**2+(V(2,1,k)-V(2,3,k))**2+&
          (V(3,1,k)-V(3,3,k))**2)
      PP4=sqrt((V(1,3,k)-V(1,4,k))**2+(V(2,3,k)-V(2,4,k))**2+&
          (V(3,3,k)-V(3,4,k))**2)
      PP5=sqrt((V(1,1,k)-V(1,4,k))**2+(V(2,1,k)-V(2,4,k))**2+&
          (V(3,1,k)-V(3,4,k))**2)
      PS1=(PP1+PP2+PP3)*0.5_DP
      PS2=(PP3+PP4+PP5)*0.5_DP
      PQ1=sqrt(PS1*(PS1-PP1)*(PS1-PP2)*(PS1-PP3))
      PQ2=sqrt(PS2*(PS2-PP3)*(PS2-PP4)*(PS2-PP5))
      PQ3=1.0_DP / (PQ1 + PQ2)
      
      ! Evalute the F2..F6 in the eight cubature points
      ! Remember: F1 = 1
      do l=1,8
        F(2,l) = CA1*P(1,l) + CB1*P(2,l) + CC1*P(3,l)
        F(3,l) = CA2*P(1,l) + CB2*P(2,l) + CC2*P(3,l)
        F(4,l) = CA3*P(1,l) + CB3*P(2,l) + CC3*P(3,l)
        F(5,l) = CD1*P(1,l)**2 + CD2*P(2,l)**2 + CD3*P(3,l)**2 &
               + CD4*P(1,l)*P(2,l) + CD5*P(1,l)*P(3,l) + CD6*P(2,l)*P(3,l)
        F(6,l) = CE1*P(1,l)**2 + CE2*P(2,l)**2 + CE3*P(3,l)**2 &
               + CE4*P(1,l)*P(2,l) + CE5*P(1,l)*P(3,l) + CE6*P(2,l)*P(3,l)
      end do

      ! Build up the matrix entries
      A(k,1)= 1.0_DP
      do l=2,6
        A(k,l) = PQ3*(PQ1*(Q8*F(l,1)+Q9*(F(l,2)+F(l,3)+F(l,4))) &
                     +PQ2*(Q8*F(l,5)+Q9*(F(l,6)+F(l,7)+F(l,8))))
      end do
               
    end do
    
    ! Now we have the coeffienct matrix - so invert it to get the
    ! coefficients of our basis polynomials.
    call mprim_invert6x6MatrixDirectDble(A,B)

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
    do k=1,6
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
    end do
    
    ! Now we can finally start to evaluate...
    dx = Dpoint(1)
    dy = Dpoint(2)
    dz = Dpoint(3)

    ! Function values
    if (BDER(DER_FUNC3D)) then
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
    end if

    ! X-derivatives
    if(BDER(DER_DERIV3D_X)) then
      Dbas(1,DER_DERIV3D_X)=2.0_DP*COB(1,1)*dx+COB(1,4)*dy+COB(1,5)*dz+COB(1,7)
      Dbas(2,DER_DERIV3D_X)=2.0_DP*COB(2,1)*dx+COB(2,4)*dy+COB(2,5)*dz+COB(2,7)
      Dbas(3,DER_DERIV3D_X)=2.0_DP*COB(3,1)*dx+COB(3,4)*dy+COB(3,5)*dz+COB(3,7)
      Dbas(4,DER_DERIV3D_X)=2.0_DP*COB(4,1)*dx+COB(4,4)*dy+COB(4,5)*dz+COB(4,7)
      Dbas(5,DER_DERIV3D_X)=2.0_DP*COB(5,1)*dx+COB(5,4)*dy+COB(5,5)*dz+COB(5,7)
      Dbas(6,DER_DERIV3D_X)=2.0_DP*COB(6,1)*dx+COB(6,4)*dy+COB(6,5)*dz+COB(6,7)
    end if

    ! Y-derivatives
    if(BDER(DER_DERIV3D_Y)) then
      Dbas(1,DER_DERIV3D_Y)=2.0_DP*COB(1,2)*dy+COB(1,4)*dx+COB(1,6)*dz+COB(1,8)
      Dbas(2,DER_DERIV3D_Y)=2.0_DP*COB(2,2)*dy+COB(2,4)*dx+COB(2,6)*dz+COB(2,8)
      Dbas(3,DER_DERIV3D_Y)=2.0_DP*COB(3,2)*dy+COB(3,4)*dx+COB(3,6)*dz+COB(3,8)
      Dbas(4,DER_DERIV3D_Y)=2.0_DP*COB(4,2)*dy+COB(4,4)*dx+COB(4,6)*dz+COB(4,8)
      Dbas(5,DER_DERIV3D_Y)=2.0_DP*COB(5,2)*dy+COB(5,4)*dx+COB(5,6)*dz+COB(5,8)
      Dbas(6,DER_DERIV3D_Y)=2.0_DP*COB(6,2)*dy+COB(6,4)*dx+COB(6,6)*dz+COB(6,8)
    end if
    
    ! Z-derivatives
    if(BDER(DER_DERIV3D_Z)) then
      Dbas(1,DER_DERIV3D_Z)=2.0_DP*COB(1,3)*dz+COB(1,5)*dx+COB(1,6)*dy+COB(1,9)
      Dbas(2,DER_DERIV3D_Z)=2.0_DP*COB(2,3)*dz+COB(2,5)*dx+COB(2,6)*dy+COB(2,9)
      Dbas(3,DER_DERIV3D_Z)=2.0_DP*COB(3,3)*dz+COB(3,5)*dx+COB(3,6)*dy+COB(3,9)
      Dbas(4,DER_DERIV3D_Z)=2.0_DP*COB(4,3)*dz+COB(4,5)*dx+COB(4,6)*dy+COB(4,9)
      Dbas(5,DER_DERIV3D_Z)=2.0_DP*COB(5,3)*dz+COB(5,5)*dx+COB(5,6)*dy+COB(5,9)
      Dbas(6,DER_DERIV3D_Z)=2.0_DP*COB(6,3)*dz+COB(6,5)*dx+COB(6,6)*dy+COB(6,9)
    end if

    ! That's it
    
  end subroutine


  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_EM30_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30_3D.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE).
  !  Dcoords(1,.)=x-coordinates,
  !  Dcoords(2,.)=y-coordinates,
  !  Dcoords(3,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i) = Coordinates of vertex i
  real(DP), dimension(:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i) = Determinant of point i
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
  ! furthermore:
  !  Dpoints(:,i) = Coordinates of point i
  real(DP), dimension(:,:), intent(IN) :: Dpoints
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
  real(DP), dimension(:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>
  real(DP), parameter :: Q1 =   0.25_DP
  real(DP), parameter :: Q2 =   1.0_DP /   3.0_DP
  real(DP), parameter :: Q3 =   0.2_DP
  real(DP), parameter :: Q4 =   0.6_DP
  real(DP), parameter :: Q8 =  -9.0_DP /  16.0_DP
  real(DP), parameter :: Q9 = 100.0_DP / 192.0_DP
  
  real(DP), dimension(3,8) :: P                       ! cubature points
  real(DP) :: PP1,PP2,PP3,PP4,PP5,PS1,PS2,PQ1,PQ2,PQ3 ! cubature weights
  real(DP), dimension(6,6) :: A, B                    ! coefficient matrices
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC1,&           ! auxiliary coefficients
              CC2,CC3,CD1,CD2,CD3,CD4,CD5,&
              CD6,CE1,CE2,CE3,CE4,CE5,CE6
  real(DP), dimension(6,10) :: COB                    ! monomial coefficients
  real(DP), dimension(3,4,6) :: V                     ! face corner vertices
  real(DP), dimension(2:6,8) :: F                     ! transformed monomials
  real(DP) :: dx,dy,dz
  integer :: i,k,l
  
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

    do k=1,6
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
      PP1=sqrt((V(1,1,k)-V(1,2,k))**2+(V(2,1,k)-V(2,2,k))**2+&
          (V(3,1,k)-V(3,2,k))**2)
      PP2=sqrt((V(1,2,k)-V(1,3,k))**2+(V(2,2,k)-V(2,3,k))**2+&
          (V(3,2,k)-V(3,3,k))**2)
      PP3=sqrt((V(1,1,k)-V(1,3,k))**2+(V(2,1,k)-V(2,3,k))**2+&
          (V(3,1,k)-V(3,3,k))**2)
      PP4=sqrt((V(1,3,k)-V(1,4,k))**2+(V(2,3,k)-V(2,4,k))**2+&
          (V(3,3,k)-V(3,4,k))**2)
      PP5=sqrt((V(1,1,k)-V(1,4,k))**2+(V(2,1,k)-V(2,4,k))**2+&
          (V(3,1,k)-V(3,4,k))**2)
      PS1=(PP1+PP2+PP3)*0.5_DP
      PS2=(PP3+PP4+PP5)*0.5_DP
      PQ1=sqrt(PS1*(PS1-PP1)*(PS1-PP2)*(PS1-PP3))
      PQ2=sqrt(PS2*(PS2-PP3)*(PS2-PP4)*(PS2-PP5))
      PQ3=1.0_DP / (PQ1 + PQ2)
      
      ! Evalute the F2..F6 in the eight cubature points
      ! Remember: F1 = 1
      do l=1,8
        F(2,l) = CA1*P(1,l) + CB1*P(2,l) + CC1*P(3,l)
        F(3,l) = CA2*P(1,l) + CB2*P(2,l) + CC2*P(3,l)
        F(4,l) = CA3*P(1,l) + CB3*P(2,l) + CC3*P(3,l)
        F(5,l) = CD1*P(1,l)**2 + CD2*P(2,l)**2 + CD3*P(3,l)**2 &
               + CD4*P(1,l)*P(2,l) + CD5*P(1,l)*P(3,l) + CD6*P(2,l)*P(3,l)
        F(6,l) = CE1*P(1,l)**2 + CE2*P(2,l)**2 + CE3*P(3,l)**2 &
               + CE4*P(1,l)*P(2,l) + CE5*P(1,l)*P(3,l) + CE6*P(2,l)*P(3,l)
      end do

      ! Build up the matrix entries
      A(k,1)= 1.0_DP
      do l=2,6
        A(k,l) = PQ3*(PQ1*(Q8*F(l,1)+Q9*(F(l,2)+F(l,3)+F(l,4))) &
                     +PQ2*(Q8*F(l,5)+Q9*(F(l,6)+F(l,7)+F(l,8))))
      end do
               
    end do
    
    ! Now we have the coeffienct matrix - so invert it
    call mprim_invert6x6MatrixDirectDble(A,B)

    ! Transform coefficiencts into monomial base
    do k=1,6
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
    end do
    
    ! Function values
    if (BDER(DER_FUNC3D)) then
      do i=1, npoints
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
      end do
    end if

    ! X-derivatives
    if(BDER(DER_DERIV3D_X)) then
      do i=1, npoints
        dx = Dpoints(1,i)
        dy = Dpoints(2,i)
        dz = Dpoints(3,i)
        Dbas(1,DER_DERIV3D_X,i)=2.0_DP*COB(1,1)*dx+COB(1,4)*dy+COB(1,5)*dz+COB(1,7)
        Dbas(2,DER_DERIV3D_X,i)=2.0_DP*COB(2,1)*dx+COB(2,4)*dy+COB(2,5)*dz+COB(2,7)
        Dbas(3,DER_DERIV3D_X,i)=2.0_DP*COB(3,1)*dx+COB(3,4)*dy+COB(3,5)*dz+COB(3,7)
        Dbas(4,DER_DERIV3D_X,i)=2.0_DP*COB(4,1)*dx+COB(4,4)*dy+COB(4,5)*dz+COB(4,7)
        Dbas(5,DER_DERIV3D_X,i)=2.0_DP*COB(5,1)*dx+COB(5,4)*dy+COB(5,5)*dz+COB(5,7)
        Dbas(6,DER_DERIV3D_X,i)=2.0_DP*COB(6,1)*dx+COB(6,4)*dy+COB(6,5)*dz+COB(6,7)
      end do
    end if

    ! Y-derivatives
    if(BDER(DER_DERIV3D_Y)) then
      do i=1, npoints
        dx = Dpoints(1,i)
        dy = Dpoints(2,i)
        dz = Dpoints(3,i)
        Dbas(1,DER_DERIV3D_Y,i)=2.0_DP*COB(1,2)*dy+COB(1,4)*dx+COB(1,6)*dz+COB(1,8)
        Dbas(2,DER_DERIV3D_Y,i)=2.0_DP*COB(2,2)*dy+COB(2,4)*dx+COB(2,6)*dz+COB(2,8)
        Dbas(3,DER_DERIV3D_Y,i)=2.0_DP*COB(3,2)*dy+COB(3,4)*dx+COB(3,6)*dz+COB(3,8)
        Dbas(4,DER_DERIV3D_Y,i)=2.0_DP*COB(4,2)*dy+COB(4,4)*dx+COB(4,6)*dz+COB(4,8)
        Dbas(5,DER_DERIV3D_Y,i)=2.0_DP*COB(5,2)*dy+COB(5,4)*dx+COB(5,6)*dz+COB(5,8)
        Dbas(6,DER_DERIV3D_Y,i)=2.0_DP*COB(6,2)*dy+COB(6,4)*dx+COB(6,6)*dz+COB(6,8)
      end do
    end if
    
    ! Z-derivatives
    if(BDER(DER_DERIV3D_Z)) then
      do i=1, npoints
        dx = Dpoints(1,i)
        dy = Dpoints(2,i)
        dz = Dpoints(3,i)
        Dbas(1,DER_DERIV3D_Z,i)=2.0_DP*COB(1,3)*dz+COB(1,5)*dx+COB(1,6)*dy+COB(1,9)
        Dbas(2,DER_DERIV3D_Z,i)=2.0_DP*COB(2,3)*dz+COB(2,5)*dx+COB(2,6)*dy+COB(2,9)
        Dbas(3,DER_DERIV3D_Z,i)=2.0_DP*COB(3,3)*dz+COB(3,5)*dx+COB(3,6)*dy+COB(3,9)
        Dbas(4,DER_DERIV3D_Z,i)=2.0_DP*COB(4,3)*dz+COB(4,5)*dx+COB(4,6)*dy+COB(4,9)
        Dbas(5,DER_DERIV3D_Z,i)=2.0_DP*COB(5,3)*dz+COB(5,5)*dx+COB(5,6)*dy+COB(5,9)
        Dbas(6,DER_DERIV3D_Z,i)=2.0_DP*COB(6,3)*dz+COB(6,5)*dx+COB(6,6)*dy+COB(6,9)
      end do
    end if
  
  ! That's it
  
  end subroutine

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_EM30_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30_3D.
  integer(I32), intent(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(IN) :: Dcoords
  
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
  real(DP), dimension(:,:,:), intent(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(IN) :: Bder
  
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
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  real(DP), dimension(:,:,:,:), intent(OUT) :: Dbas
!</output>

! </subroutine>
  real(DP), parameter :: Q1 =   0.25_DP
  real(DP), parameter :: Q2 =   1.0_DP /   3.0_DP
  real(DP), parameter :: Q3 =   0.2_DP
  real(DP), parameter :: Q4 =   0.6_DP
  real(DP), parameter :: Q8 =  -9.0_DP /  16.0_DP
  real(DP), parameter :: Q9 = 100.0_DP / 192.0_DP
  
  real(DP), dimension(3,8) :: P                       ! cubature points
  real(DP) :: PP1,PP2,PP3,PP4,PP5,PS1,PS2,PQ1,PQ2,PQ3 ! cubature weights
  real(DP), dimension(6,6) :: A, B                    ! coefficient matrices
  real(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC1,&           ! auxiliary coefficients
              CC2,CC3,CD1,CD2,CD3,CD4,CD5,&
              CD6,CE1,CE2,CE3,CE4,CE5,CE6
  real(DP), dimension(6,10) :: COB                    ! monomial coefficients
  real(DP), dimension(3,4,6) :: V                     ! face corner vertices
  real(DP), dimension(2:6,8) :: F                     ! transformed monomials
  real(DP) :: dx,dy,dz
  integer :: i,j,k,l
  
  ! Loop over all elements
  do j=1, nelements
  
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

    do k=1,6
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
      PP1=sqrt((V(1,1,k)-V(1,2,k))**2+(V(2,1,k)-V(2,2,k))**2+&
          (V(3,1,k)-V(3,2,k))**2)
      PP2=sqrt((V(1,2,k)-V(1,3,k))**2+(V(2,2,k)-V(2,3,k))**2+&
          (V(3,2,k)-V(3,3,k))**2)
      PP3=sqrt((V(1,1,k)-V(1,3,k))**2+(V(2,1,k)-V(2,3,k))**2+&
          (V(3,1,k)-V(3,3,k))**2)
      PP4=sqrt((V(1,3,k)-V(1,4,k))**2+(V(2,3,k)-V(2,4,k))**2+&
          (V(3,3,k)-V(3,4,k))**2)
      PP5=sqrt((V(1,1,k)-V(1,4,k))**2+(V(2,1,k)-V(2,4,k))**2+&
          (V(3,1,k)-V(3,4,k))**2)
      PS1=(PP1+PP2+PP3)*0.5_DP
      PS2=(PP3+PP4+PP5)*0.5_DP
      PQ1=sqrt(PS1*(PS1-PP1)*(PS1-PP2)*(PS1-PP3))
      PQ2=sqrt(PS2*(PS2-PP3)*(PS2-PP4)*(PS2-PP5))
      PQ3=1.0_DP / (PQ1 + PQ2)
      
      ! Evalute the F2..F6 in the eight cubature points
      ! Remember: F1 = 1
      do l=1,8
        F(2,l) = CA1*P(1,l) + CB1*P(2,l) + CC1*P(3,l)
        F(3,l) = CA2*P(1,l) + CB2*P(2,l) + CC2*P(3,l)
        F(4,l) = CA3*P(1,l) + CB3*P(2,l) + CC3*P(3,l)
        F(5,l) = CD1*P(1,l)**2 + CD2*P(2,l)**2 + CD3*P(3,l)**2 &
               + CD4*P(1,l)*P(2,l) + CD5*P(1,l)*P(3,l) + CD6*P(2,l)*P(3,l)
        F(6,l) = CE1*P(1,l)**2 + CE2*P(2,l)**2 + CE3*P(3,l)**2 &
               + CE4*P(1,l)*P(2,l) + CE5*P(1,l)*P(3,l) + CE6*P(2,l)*P(3,l)
      end do

      ! Build up the matrix entries
      A(k,1)= 1.0_DP
      do l=2,6
        A(k,l) = PQ3*(PQ1*(Q8*F(l,1)+Q9*(F(l,2)+F(l,3)+F(l,4))) &
                     +PQ2*(Q8*F(l,5)+Q9*(F(l,6)+F(l,7)+F(l,8))))
      end do
               
    end do
    
    ! Now we have the coeffienct matrix - so invert it
    call mprim_invert6x6MatrixDirectDble(A,B)

    ! Transform coefficiencts into monomial base
    do k=1,6
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
    end do
    
    ! Function values
    if (BDER(DER_FUNC3D)) then
      do i=1, npoints
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
      end do
    end if

    ! X-derivatives
    if(BDER(DER_DERIV3D_X)) then
      do i=1, npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        Dbas(1,DER_DERIV3D_X,i,j)=2.0_DP*COB(1,1)*dx+COB(1,4)*dy+COB(1,5)*dz+COB(1,7)
        Dbas(2,DER_DERIV3D_X,i,j)=2.0_DP*COB(2,1)*dx+COB(2,4)*dy+COB(2,5)*dz+COB(2,7)
        Dbas(3,DER_DERIV3D_X,i,j)=2.0_DP*COB(3,1)*dx+COB(3,4)*dy+COB(3,5)*dz+COB(3,7)
        Dbas(4,DER_DERIV3D_X,i,j)=2.0_DP*COB(4,1)*dx+COB(4,4)*dy+COB(4,5)*dz+COB(4,7)
        Dbas(5,DER_DERIV3D_X,i,j)=2.0_DP*COB(5,1)*dx+COB(5,4)*dy+COB(5,5)*dz+COB(5,7)
        Dbas(6,DER_DERIV3D_X,i,j)=2.0_DP*COB(6,1)*dx+COB(6,4)*dy+COB(6,5)*dz+COB(6,7)
      end do
    end if

    ! Y-derivatives
    if(BDER(DER_DERIV3D_Y)) then
      do i=1, npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        Dbas(1,DER_DERIV3D_Y,i,j)=2.0_DP*COB(1,2)*dy+COB(1,4)*dx+COB(1,6)*dz+COB(1,8)
        Dbas(2,DER_DERIV3D_Y,i,j)=2.0_DP*COB(2,2)*dy+COB(2,4)*dx+COB(2,6)*dz+COB(2,8)
        Dbas(3,DER_DERIV3D_Y,i,j)=2.0_DP*COB(3,2)*dy+COB(3,4)*dx+COB(3,6)*dz+COB(3,8)
        Dbas(4,DER_DERIV3D_Y,i,j)=2.0_DP*COB(4,2)*dy+COB(4,4)*dx+COB(4,6)*dz+COB(4,8)
        Dbas(5,DER_DERIV3D_Y,i,j)=2.0_DP*COB(5,2)*dy+COB(5,4)*dx+COB(5,6)*dz+COB(5,8)
        Dbas(6,DER_DERIV3D_Y,i,j)=2.0_DP*COB(6,2)*dy+COB(6,4)*dx+COB(6,6)*dz+COB(6,8)
      end do
    end if
    
    ! Z-derivatives
    if(BDER(DER_DERIV3D_Z)) then
      do i=1, npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        Dbas(1,DER_DERIV3D_Z,i,j)=2.0_DP*COB(1,3)*dz+COB(1,5)*dx+COB(1,6)*dy+COB(1,9)
        Dbas(2,DER_DERIV3D_Z,i,j)=2.0_DP*COB(2,3)*dz+COB(2,5)*dx+COB(2,6)*dy+COB(2,9)
        Dbas(3,DER_DERIV3D_Z,i,j)=2.0_DP*COB(3,3)*dz+COB(3,5)*dx+COB(3,6)*dy+COB(3,9)
        Dbas(4,DER_DERIV3D_Z,i,j)=2.0_DP*COB(4,3)*dz+COB(4,5)*dx+COB(4,6)*dy+COB(4,9)
        Dbas(5,DER_DERIV3D_Z,i,j)=2.0_DP*COB(5,3)*dz+COB(5,5)*dx+COB(5,6)*dy+COB(5,9)
        Dbas(6,DER_DERIV3D_Z,i,j)=2.0_DP*COB(6,3)*dz+COB(6,5)*dx+COB(6,6)*dy+COB(6,9)
      end do
    end if

  end do ! element loop
  
  ! That's it
  end subroutine

end module