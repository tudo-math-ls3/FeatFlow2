!##############################################################################
!# ****************************************************************************
!# <name> element_prism3d </name>
!# ****************************************************************************
!# 
!# <purpose>
!# This module contains the implementations of the 3D prism basis
!# functions.
!#
!# </purpose>
!##############################################################################

module element_prism3d

  use basicgeometry
  use derivatives
  use elementbase
  use fsystem
  use perfconfig

  implicit none
  
  private
  
  public :: elem_R0_3D 
  public :: elem_R0_3D_mult 
  public :: elem_R0_3D_sim 
  public :: elem_R1_3D 
  public :: elem_R1_3D_mult 
  public :: elem_R1_3D_sim 

contains

!**************************************************************************
! Element subroutines for parametric 3D R0 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  pure subroutine elem_R0_3D (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_R0_3D.
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
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  real(DP), intent(in) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  ! Dpoint(3) = z-coordinate
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
    
  ! R0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC) = 1.0_DP

  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_R0_3D_mult (celement, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_R0_3D.
  integer(I32), intent(in)  :: celement
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
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
  !  Djac(9,i) = J_I(3,3)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:), intent(in) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates,
  ! Dpoints(3,.)=z-coordinates.
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
  
  ! R0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:) = 1.0_DP

  end subroutine 

  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_R0_3D_sim (celement, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_R0_3D.
  integer(I32), intent(in)  :: celement
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
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
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:,:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
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
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
  
  !</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! R0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  end subroutine 

!**************************************************************************
! Element subroutines for parametric R1 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  pure subroutine elem_R1_3D (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_R1_3D.
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
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), intent(in) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate,
  ! Dpoint(3) = z-coordinate
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

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,NDIM3D) :: Dhelp
  real(DP) :: dx,dy,dz, djx, djy, djz

  real(DP) :: dxj !auxiliary variable
  
  ! The R1 element is specified by six polynomials on the reference element.
  ! These six polynomials are:
  !
  !  P1(x,y,z) = 1/2 (1-x-y) (1-z)
  !  P2(x,y,z) = 1/2 x (1-z)
  !  P3(x,y,z) = 1/2 y (1-z)
  !  P4(x,y,z) = 1/2 (1-x-y) (1+z)
  !  P5(x,y,z) = 1/2 x (1+z)
  !  P6(x,y,z) = 1/2 y (1+z)
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,...,X6 the six corners of the reference prism.
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
  dz = Dpoint(3)
    
  ! Remark: The R1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC3D)) then
    Dbas(1,DER_FUNC3D) = 0.5_DP*(1.0_DP-dx-dy)*(1.0_DP-dz)
    Dbas(2,DER_FUNC3D) = 0.5_DP*dx*(1.0_DP-dz)
    Dbas(3,DER_FUNC3D) = 0.5_DP*dy*(1.0_DP-dz)
    Dbas(4,DER_FUNC3D) = 0.5_DP*(1.0_DP-dx-dy)*(1.0_DP+dz)
    Dbas(5,DER_FUNC3D) = 0.5_DP*dx*(1.0_DP+dz)
    Dbas(6,DER_FUNC3D) = 0.5_DP*dy*(1.0_DP+dz)
!  endif
  
  ! If x-, y- or z-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
!      (Bder(DER_DERIV3D_Z))) then
    dxj = 0.5_DP / ddetj
    
    ! x-, y- and z-derivatives on reference element
    Dhelp(1,1) =-(1.0_DP-dz)
    Dhelp(2,1) = (1.0_DP-dz)
    Dhelp(3,1) = 0.0_DP
    Dhelp(4,1) =-(1.0_DP+dz)
    Dhelp(5,1) = (1.0_DP+dz)
    Dhelp(6,1) = 0.0_DP
    Dhelp(1,2) =-(1.0_DP-dz)
    Dhelp(2,2) = 0.0_DP
    Dhelp(3,2) = (1.0_DP-dz)
    Dhelp(4,2) =-(1.0_DP+dz)
    Dhelp(5,2) = 0.0_DP
    Dhelp(6,2) = (1.0_DP+dz)
    Dhelp(1,3) =-(1.0_DP-dx-dy)
    Dhelp(2,3) =-dx
    Dhelp(3,3) =-dy
    Dhelp(4,3) = (1.0_DP-dx-dy)
    Dhelp(5,3) = dx
    Dhelp(6,3) = dy
      
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
!    endif
!  endif
    
  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  pure subroutine elem_R1_3D_mult (celement, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_R1_3D.
  integer(I32), intent(in)  :: celement
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
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
  !  Djac(9,i) = J_I(3,3)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  real(DP), dimension(:,:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:), intent(in) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
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
  real(DP), dimension(6,NDIM3D) :: Dhelp

  real(DP) :: djx, djy, djz,dxj
  integer :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The R1-element always computes function value and 1st derivatives.
  ! That is even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC3D,i) = 0.5_DP*(1.0_DP-Dpoints(1,i)-Dpoints(2,i))&
                                   *(1.0_DP-Dpoints(3,i))
      Dbas(2,DER_FUNC3D,i) = 0.5_DP*Dpoints(1,i)*(1.0_DP-Dpoints(3,i))
      Dbas(3,DER_FUNC3D,i) = 0.5_DP*Dpoints(2,i)*(1.0_DP-Dpoints(3,i))
      Dbas(5,DER_FUNC3D,i) = 0.5_DP*(1.0_DP-Dpoints(1,i)-Dpoints(2,i))&
                                   *(1.0_DP+Dpoints(3,i))
      Dbas(6,DER_FUNC3D,i) = 0.5_DP*Dpoints(1,i)*(1.0_DP+Dpoints(3,i))
      Dbas(7,DER_FUNC3D,i) = 0.5_DP*Dpoints(2,i)*(1.0_DP+Dpoints(3,i))

    end do
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR.&
!      (Bder(DER_DERIV3D_Z))) THEN
    
    !x-, y- and z-derivatives on reference element
    do i = 1, npoints
      Dhelp(1,1) =-(1.0_DP-Dpoints(3,i))
      Dhelp(2,1) = (1.0_DP-Dpoints(3,i))
      Dhelp(3,1) = 0.0_DP
      Dhelp(4,1) =-(1.0_DP+Dpoints(3,i))
      Dhelp(5,1) = (1.0_DP+Dpoints(3,i))
      Dhelp(6,1) = 0.0_DP
      Dhelp(1,2) =-(1.0_DP-Dpoints(3,i))
      Dhelp(2,2) = 0.0_DP
      Dhelp(3,2) = (1.0_DP-Dpoints(3,i))
      Dhelp(4,2) =-(1.0_DP+Dpoints(3,i))
      Dhelp(5,2) = 0.0_DP
      Dhelp(6,2) = (1.0_DP+Dpoints(3,i))
      Dhelp(1,3) =-(1.0_DP-Dpoints(1,i)-Dpoints(2,i))
      Dhelp(2,3) =-Dpoints(1,i)
      Dhelp(3,3) =-Dpoints(2,i)
      Dhelp(4,3) = (1.0_DP-Dpoints(1,i)-Dpoints(2,i))
      Dhelp(5,3) = Dpoints(1,i)
      Dhelp(6,3) = Dpoints(2,i)
      
      dxj = 0.5_DP / Ddetj(i)
      
      ! x-derivatives on current element
      djx = Djac(5,i)*Djac(9,i) - Djac(6,i)*Djac(8,i)
      djy = Djac(8,i)*Djac(3,i) - Djac(2,i)*Djac(9,i)
      djz = Djac(2,i)*Djac(6,i) - Djac(5,i)*Djac(3,i)
      Dbas(1,DER_DERIV3D_X,i) = Dxj * &
          (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_X,i) = Dxj * &
          (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_X,i) = Dxj * &
          (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_X,i) = Dxj * &
          (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_X,i) = Dxj * &
          (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_X,i) = Dxj * &
          (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
    
      ! y-derivatives on current element
      djx = Djac(7,i)*Djac(6,i) - Djac(4,i)*Djac(9,i)
      djy = Djac(1,i)*Djac(9,i) - Djac(7,i)*Djac(3,i)
      djz = Djac(4,i)*Djac(3,i) - Djac(1,i)*Djac(6,i)
      Dbas(1,DER_DERIV3D_Y,i) = Dxj * &
          (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Y,i) = Dxj * &
          (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Y,i) = Dxj * &
          (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Y,i) = Dxj * &
          (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Y,i) = Dxj * &
          (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Y,i) = Dxj * &
          (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))

      ! z-derivatives on current element
      djx = Djac(4,i)*Djac(8,i) - Djac(7,i)*Djac(5,i)
      djy = Djac(7,i)*Djac(2,i) - Djac(1,i)*Djac(8,i)
      djz = Djac(1,i)*Djac(5,i) - Djac(4,i)*Djac(2,i)
      Dbas(1,DER_DERIV3D_Z,i) = Dxj * &
          (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Z,i) = Dxj * &
          (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Z,i) = Dxj * &
          (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Z,i) = Dxj * &
          (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Z,i) = Dxj * &
          (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Z,i) = Dxj * &
          (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
    end do
!  ENDIF
    
  end subroutine 

  !************************************************************************
  
!<subroutine>  

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_R1_3D_sim (celement, Dcoords, Djac, Ddetj, &
                             Bder, Dbas, npoints, nelements, &
                             Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_R1_3D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
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
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
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
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(out) :: Dbas
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,NDIM3D) :: Dhelp
  real(DP) :: djx, djy, djz, dxj
  
  integer :: i   ! point counter
  integer :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  if (Bder(DER_FUNC3D)) then
  
    !$omp parallel do default(shared) private(i) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
    
      do i=1,npoints
        Dbas(1,DER_FUNC3D,i,j) = 0.25_DP*(1.0_DP-Dpoints(1,i,j)-Dpoints(2,i,j))&
                                        *(1.0_DP-Dpoints(3,i,j))
        Dbas(2,DER_FUNC3D,i,j) = 0.25_DP*Dpoints(1,i,j)*(1.0_DP-Dpoints(3,i,j))
        Dbas(3,DER_FUNC3D,i,j) = 0.25_DP*Dpoints(2,i,j)*(1.0_DP-Dpoints(3,i,j))
        Dbas(4,DER_FUNC3D,i,j) = 0.25_DP*(1.0_DP-Dpoints(1,i,j)-Dpoints(2,i,j))&
                                        *(1.0_DP+Dpoints(3,i,j))
        Dbas(5,DER_FUNC3D,i,j) = 0.25_DP*Dpoints(1,i,j)*(1.0_DP+Dpoints(3,i,j))
        Dbas(6,DER_FUNC3D,i,j) = 0.25_DP*Dpoints(2,i,j)*(1.0_DP+Dpoints(3,i,j))
      end do
      
    end do
    !$omp end parallel do
    
  end if
    
  !if x-, y- or z-derivatives are desired
  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
      (Bder(DER_DERIV3D_Z))) then
  
    !$omp parallel do default(shared) private(i,Dxj,Dhelp,djx,djy,djz) &
    !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
    do j=1,nelements
      
      ! x-, y- and z-derivatives on reference element
      do i=1,npoints
        Dhelp(1,1) =-(1.0_DP-Dpoints(3,i,j))
        Dhelp(2,1) = (1.0_DP-Dpoints(3,i,j))
        Dhelp(3,1) = 0.0_DP
        Dhelp(4,1) =-(1.0_DP+Dpoints(3,i,j))
        Dhelp(5,1) = (1.0_DP+Dpoints(3,i,j))
        Dhelp(6,1) = 0.0_DP
        Dhelp(1,2) =-(1.0_DP-Dpoints(3,i,j))
        Dhelp(2,2) = 0.0_DP
        Dhelp(3,2) = (1.0_DP-Dpoints(3,i,j))
        Dhelp(4,2) =-(1.0_DP+Dpoints(3,i,j))
        Dhelp(5,2) = 0.0_DP
        Dhelp(6,2) = (1.0_DP+Dpoints(3,i,j))
        Dhelp(1,3) =-(1.0_DP-Dpoints(1,i,j)-Dpoints(2,i,j))
        Dhelp(2,3) =-Dpoints(1,i,j)
        Dhelp(3,3) =-Dpoints(2,i,j)
        Dhelp(4,3) =+(1.0_DP-Dpoints(1,i,j)-Dpoints(2,i,j))
        Dhelp(5,3) =+Dpoints(1,i,j)
        Dhelp(6,3) =+Dpoints(2,i,j)

        dxj = 0.25_DP / Ddetj(i,j)
        
        ! x-derivatives on current element
        djx = Djac(5,i,j)*Djac(9,i,j) - Djac(6,i,j)*Djac(8,i,j)
        djy = Djac(8,i,j)*Djac(3,i,j) - Djac(2,i,j)*Djac(9,i,j)
        djz = Djac(2,i,j)*Djac(6,i,j) - Djac(5,i,j)*Djac(3,i,j)
        Dbas(1,DER_DERIV3D_X,i,j) = Dxj * &
            (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
        Dbas(2,DER_DERIV3D_X,i,j) = Dxj * &
            (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
        Dbas(3,DER_DERIV3D_X,i,j) = Dxj * &
            (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
        Dbas(4,DER_DERIV3D_X,i,j) = Dxj * &
            (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
        Dbas(5,DER_DERIV3D_X,i,j) = Dxj * &
            (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
        Dbas(6,DER_DERIV3D_X,i,j) = Dxj * &
            (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
      
        ! y-derivatives on current element
        djx = Djac(7,i,j)*Djac(6,i,j) - Djac(4,i,j)*Djac(9,i,j)
        djy = Djac(1,i,j)*Djac(9,i,j) - Djac(7,i,j)*Djac(3,i,j)
        djz = Djac(4,i,j)*Djac(3,i,j) - Djac(1,i,j)*Djac(6,i,j)
        Dbas(1,DER_DERIV3D_Y,i,j) = Dxj * &
            (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
        Dbas(2,DER_DERIV3D_Y,i,j) = Dxj * &
            (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
        Dbas(3,DER_DERIV3D_Y,i,j) = Dxj * &
            (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
        Dbas(4,DER_DERIV3D_Y,i,j) = Dxj * &
            (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
        Dbas(5,DER_DERIV3D_Y,i,j) = Dxj * &
            (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
        Dbas(6,DER_DERIV3D_Y,i,j) = Dxj * &
            (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))

        ! z-derivatives on current element
        djx = Djac(4,i,j)*Djac(8,i,j) - Djac(7,i,j)*Djac(5,i,j)
        djy = Djac(7,i,j)*Djac(2,i,j) - Djac(1,i,j)*Djac(8,i,j)
        djz = Djac(1,i,j)*Djac(5,i,j) - Djac(4,i,j)*Djac(2,i,j)
        Dbas(1,DER_DERIV3D_Z,i,j) = Dxj * &
            (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
        Dbas(2,DER_DERIV3D_Z,i,j) = Dxj * &
            (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
        Dbas(3,DER_DERIV3D_Z,i,j) = Dxj * &
            (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
        Dbas(4,DER_DERIV3D_Z,i,j) = Dxj * &
            (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
        Dbas(5,DER_DERIV3D_Z,i,j) = Dxj * &
            (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
        Dbas(6,DER_DERIV3D_Z,i,j) = Dxj * &
            (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))

      end do
    end do
    !$omp end parallel do
      
  end if
    
  end subroutine 

end module
