!##############################################################################
!# ****************************************************************************
!# <name> element_line1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the implementations of the 1D line basis functions.
!#
!# </purpose>
!##############################################################################

module element_line1d

!$ use omp_lib
  use fsystem
  use elementbase
  use derivatives
  use perfconfig

  implicit none

  private

  public :: elem_P0_1D
  public :: elem_P0_1D_mult
  public :: elem_P0_1D_sim
  public :: elem_P1_1D
  public :: elem_P1_1D_mult
  public :: elem_P1_1D_sim
  public :: elem_P2_1D
  public :: elem_P2_1D_mult
  public :: elem_P2_1D_sim
  public :: elem_S31_1D
  public :: elem_S31_1D_mult
  public :: elem_S31_1D_sim
  public :: elem_eval_P1_1D
  public :: elem_eval_P2_1D
  public :: elem_eval_S31_1D
  public :: elem_eval_PN_1D
  public :: elem_DG_T0_1D
  public :: elem_DG_T0_1D_mult
  public :: elem_DG_T0_1D_sim
  public :: elem_DG_T1_1D
  public :: elem_DG_T1_1D_mult
  public :: elem_DG_T1_1D_sim
  public :: elem_DG_T2_1D
  public :: elem_DG_T2_1D_mult
  public :: elem_DG_T2_1D_sim

contains

!**************************************************************************
! Element subroutines for parametric P0_1D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_P0_1D (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P0_1D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
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
  real(DP), dimension(1), intent(in) :: Dpoint
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

  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC1D) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P0_1D_mult (celement, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P0_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
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

  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC1D,:) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P0_1D_sim (celement, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P0_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: ddetj

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

  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC1D,:,:) = 1.0_DP

  end subroutine

!**************************************************************************
! Element subroutines for parametric P1_1D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_P1_1D (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P1_1D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
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
  real(DP), dimension(1), intent(in) :: Dpoint
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
  ! That is even faster than when using two IF commands for preventing
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

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P1_1D_mult (celement, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
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

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1_1D-element always computes function value and 1st derivatives
  ! That is even faster than when using two IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC1D)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC1D,i) = 0.5_DP * (1_DP - Dpoints(1,i))
      Dbas(2,DER_FUNC1D,i) = 0.5_DP * (1_DP + Dpoints(1,i))
    !end do
  !ENDIF

  !if x-derivatives are desired
!  IF (Bder(DER_DERIV1D_X)) THEN
    !do i=1,npoints
      Dbas(1,DER_DERIV1D_X,i) = -0.5_DP / Ddetj(i)
      Dbas(2,DER_DERIV1D_X,i) = 0.5_DP / Ddetj(i)
    end do
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_P1_1D_sim (celement, Dcoords, Djac, Ddetj, &
                             Bder, Dbas, npoints, nelements, &
                             Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
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

  integer :: i   ! point counter
  integer :: j   ! element counter

    !if function values are desired
    if (Bder(DER_FUNC1D)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_FUNC1D,i,j) = 0.5_DP * (1.0_DP - Dpoints(1,i,j))
          Dbas(2,DER_FUNC1D,i,j) = 0.5_DP * (1.0_DP + Dpoints(1,i,j))
        end do

      end do
      !$omp end parallel do

    end if

    !if x-derivatives are desired
    if (Bder(DER_DERIV1D_X)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_DERIV1D_X,i,j) = -0.5 / Ddetj(i,j)
          Dbas(2,DER_DERIV1D_X,i,j) = 0.5 / Ddetj(i,j)
        end do

      end do
      !$omp end parallel do

    end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric P2_1D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_P2_1D (celement, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P2_1D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
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
  real(DP), dimension(1), intent(in) :: Dpoint
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

  real(DP) :: d

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
  ! That is even faster than when using two IF commands for preventing
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

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_P2_1D_mult (celement, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P2_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
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

  real(DP) :: d
  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P2_1D-element always computes function value and 1st derivatives
  ! That is even faster than when using two IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC1D)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC1D,i) = 0.5_DP * Dpoints(1,i) * (Dpoints(1,i) - 1.0_DP)
      Dbas(2,DER_FUNC1D,i) = 0.5_DP * Dpoints(1,i) * (Dpoints(1,i) + 1.0_DP)
      Dbas(3,DER_FUNC1D,i) = 1.0_DP - Dpoints(1,i)**2
    !end do
  !ENDIF

  !if x-derivatives are desired
!  IF (Bder(DER_DERIV1D_X)) THEN
    !do i=1,npoints
      d = 1.0_DP / Ddetj(i)
      Dbas(1,DER_DERIV1D_X,i) = (Dpoints(1,i) - 0.5_DP) * d
      Dbas(2,DER_DERIV1D_X,i) = (Dpoints(1,i) + 0.5_DP) * d
      Dbas(3,DER_DERIV1D_X,i) = -2.0_DP * Dpoints(1,i) * d
    end do
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_P2_1D_sim (celement, Dcoords, Djac, Ddetj, &
                             Bder, Dbas, npoints, nelements, &
                             Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P2_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
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

  real(DP) :: d
  integer :: i   ! point counter
  integer :: j   ! element counter

    !if function values are desired
    if (Bder(DER_FUNC1D)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_FUNC1D,i,j) = 0.5_DP*Dpoints(1,i,j)*(Dpoints(1,i,j) - 1.0_DP)
          Dbas(2,DER_FUNC1D,i,j) = 0.5_DP*Dpoints(1,i,j)*(Dpoints(1,i,j) + 1.0_DP)
          Dbas(3,DER_FUNC1D,i,j) = 1.0_DP - Dpoints(1,i,j)**2
        end do

      end do
      !$omp end parallel do

    end if

    !if x-derivatives are desired
    if (Bder(DER_DERIV1D_X)) then

      !$omp parallel do default(shared) private(i,d) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          d = 1.0_DP / Ddetj(i,j)
          Dbas(1,DER_DERIV1D_X,i,j) = (Dpoints(1,i,j) - 0.5_DP) * d
          Dbas(2,DER_DERIV1D_X,i,j) = (Dpoints(1,i,j) + 0.5_DP) * d
          Dbas(3,DER_DERIV1D_X,i,j) = -2.0_DP * Dpoints(1,i,j) * d
        end do

      end do
      !$omp end parallel do

    end if

  end subroutine

!**************************************************************************
! Element subroutines for parametric S31_1D element.
! The routines are defines with the F95 PURE statement as they work
! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_S31_1D (celement, Dcoords, Djac, ddetj, Bder, &
                               Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_S31_1D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
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
  real(DP), dimension(1), intent(in) :: Dpoint
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

  ! auxiliary variable
  real(DP) :: dxj

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
    Dbas(3,DER_DERIV1D_X) = 0.25_DP*(3.0_DP*Dpoint(1) + 1.0_DP)*(Dpoint(1) - 1.0_DP)
    Dbas(4,DER_DERIV1D_X) = 0.25_DP*(3.0_DP*Dpoint(1) - 1.0_DP)*(Dpoint(1) + 1.0_DP)
!  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_S31_1D_mult (celement, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_S31_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
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

  integer :: i   ! point counter
  real(DP) :: dxj

    !if function values are desired
    if (Bder(DER_FUNC1D)) then
      do i=1,npoints
        Dbas(1,DER_FUNC1D,i) = 0.25_DP*Dpoints(1,i)*&
                              (Dpoints(1,i)**2 - 3.0_DP) + 0.5_DP
        Dbas(2,DER_FUNC1D,i) = 0.25_DP*Dpoints(1,i)*&
                              (3.0_DP - Dpoints(1,i)**2) + 0.5_DP
        Dbas(3,DER_FUNC1D,i) = Ddetj(i)*0.25*(Dpoints(1,i) + 1.0_DP)*&
                              (Dpoints(1,i) - 1.0_DP)**2
        Dbas(4,DER_FUNC1D,i) = Ddetj(i)*0.25*(Dpoints(1,i) - 1.0_DP)*&
                              (Dpoints(1,i) + 1.0_DP)**2
      end do
    endif

    !if x-derivatives are desired
    if (Bder(DER_DERIV1D_X)) then
      do i=1,npoints
        dxj = 1.0_DP / Ddetj(i)
        Dbas(1,DER_DERIV1D_X,i) = 0.75_DP*(Dpoints(1,i)**2 - 1.0_DP)*dxj
        Dbas(2,DER_DERIV1D_X,i) = 0.75_DP*(1.0_DP - Dpoints(1,i)**2)*dxj
        Dbas(3,DER_DERIV1D_X,i) = 0.25_DP*(3.0_DP*Dpoints(1,i) + 1.0_DP)*(Dpoints(1,i) - 1.0_DP)
        Dbas(4,DER_DERIV1D_X,i) = 0.25_DP*(3.0_DP*Dpoints(1,i) - 1.0_DP)*(Dpoints(1,i) + 1.0_DP)
      end do
    endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_S31_1D_sim (celement, Dcoords, Djac, Ddetj, &
                              Bder, Dbas, npoints, nelements, &
                              Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
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

  integer :: i   ! point counter
  integer :: j   ! element counter
  real(DP) :: dxj

    !if function values are desired
    if (Bder(DER_FUNC1D)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_FUNC1D,i,j) = 0.25_DP*Dpoints(1,i,j)*&
                                  (Dpoints(1,i,j)**2 - 3.0_DP) + 0.5_DP
          Dbas(2,DER_FUNC1D,i,j) = 0.25_DP*Dpoints(1,i,j)*&
                                  (3.0_DP - Dpoints(1,i,j)**2) + 0.5_DP
          Dbas(3,DER_FUNC1D,i,j) = Ddetj(i,j)*0.25*(Dpoints(1,i,j) + 1.0_DP)*&
                                  (Dpoints(1,i,j) - 1.0_DP)**2
          Dbas(4,DER_FUNC1D,i,j) = Ddetj(i,j)*0.25*(Dpoints(1,i,j) - 1.0_DP)*&
                                  (Dpoints(1,i,j) + 1.0_DP)**2
        end do

      end do
      !$omp end parallel do

    end if

    !if x-derivatives are desired
    if (Bder(DER_DERIV1D_X)) then

      !$omp parallel do default(shared) private(i,dxj) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          dxj = 1.0_DP / Ddetj(i,j)
          Dbas(1,DER_DERIV1D_X,i,j) = 0.75_DP*(Dpoints(1,i,j)**2 - 1.0_DP)*dxj
          Dbas(2,DER_DERIV1D_X,i,j) = 0.75_DP*(1.0_DP - Dpoints(1,i,j)**2)*dxj
          Dbas(3,DER_DERIV1D_X,i,j) = 0.25_DP*(3.0_DP*Dpoints(1,i,j) + 1.0_DP)*(Dpoints(1,i,j) - 1.0_DP)
          Dbas(4,DER_DERIV1D_X,i,j) = 0.25_DP*(3.0_DP*Dpoints(1,i,j) - 1.0_DP)*(Dpoints(1,i,j) + 1.0_DP)
        end do

      end do
      !$omp end parallel do

    end if

  end subroutine

  !*****************************************************************************
  ! Element subroutines for parametric DG_T0_1D element.
  ! The routines are defines with the F95 PURE statement as they work
  ! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_DG_T0_1D (celement, Dcoords, Djac, ddetj, &
                                 Bder, Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T0_1D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
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
  real(DP), dimension(1), intent(in) :: Dpoint
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

  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC1D) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T0_1D_mult (celement, Dcoords, Djac, Ddetj, &
                                      Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T0_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
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

  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC1D,:) = 1.0_DP

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T0_1D_sim (celement, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T0_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:,:), intent(in) :: Djac

  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  real(DP), dimension(:,:), intent(in) :: ddetj

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

  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC1D,:,:) = 1.0_DP

  end subroutine

  !*****************************************************************************
  ! Element subroutines for parametric DG_T1_1D element.
  ! The routines are defines with the F95 PURE statement as they work
  ! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_DG_T1_1D (celement, Dcoords, Djac, ddetj, Bder, &
                                 Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T1_1D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
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
  real(DP), dimension(1), intent(in) :: Dpoint
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

  ! The DG_T1_1D element is specified by two polynomials on the reference element.
  ! These two polynomials are:
  !
  !  P1(x) = 1
  !  P2(x) = x
  !
  ! The first one is obviously the constant polynomial and the second one
  ! has vanishing mean-value on the reference element [-1,1].

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T1_1D-element always computes function value and 1st derivative.
  ! That is even faster than when using two IF commands for preventing
  ! the computation of one of the values!

  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC1D)) then
    Dbas(1,DER_FUNC1D) = 1.0_DP
    Dbas(2,DER_FUNC1D) = Dpoint(1)
!  endif

  ! If x-derivatives are desired, calculate them.
  ! Since DG_T1_1D is linear, the first derivative is constant!
!  if (Bder(DER_DERIV1D_X)) then
    Dbas(1,DER_DERIV1D_X) = 0.0_DP
    Dbas(2,DER_DERIV1D_X) = 1.0_DP / ddetj
!  endif

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T1_1D_mult (celement, Dcoords, Djac, Ddetj, &
                                      Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T1_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
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

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T1_1D-element always computes function value and 1st derivatives
  ! That is even faster than when using two IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC1D)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC1D,i) = 1.0_DP
      Dbas(2,DER_FUNC1D,i) = Dpoints(1,i)
    !end do
  !ENDIF

  !if x-derivatives are desired
!  IF (Bder(DER_DERIV1D_X)) THEN
    !do i=1,npoints
      Dbas(1,DER_DERIV1D_X,i) = 0.0_DP
      Dbas(2,DER_DERIV1D_X,i) = 1.0_DP / Ddetj(i)
    end do
!  ENDIF

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_DG_T1_1D_sim (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, nelements, &
                                Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T1_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
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

  integer :: i   ! point counter
  integer :: j   ! element counter

    !if function values are desired
    if (Bder(DER_FUNC1D)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_FUNC1D,i,j) = 1.0_DP
          Dbas(2,DER_FUNC1D,i,j) = Dpoints(1,i,j)
        end do

      end do
      !$omp end parallel do

    end if

    !if x-derivatives are desired
    if (Bder(DER_DERIV1D_X)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_DERIV1D_X,i,j) = 0.0_DP
          Dbas(2,DER_DERIV1D_X,i,j) = 1.0_DP / Ddetj(i,j)
        end do

      end do
      !$omp end parallel do

    end if

  end subroutine

  !*****************************************************************************
  ! Element subroutines for parametric DG_T2_1D element.
  ! The routines are defines with the F95 PURE statement as they work
  ! only on the parameters; helps some compilers in optimisation.

!<subroutine>

  pure subroutine elem_DG_T2_1D (celement, Dcoords, Djac, ddetj, Bder, &
                                 Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_DG_T2_1D.
  integer(I32), intent(in)  :: celement

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
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
  real(DP), dimension(1), intent(in) :: Dpoint
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

  ! The DG_T2_1D element is specified by three polynomials on the reference element.
  ! These three polynomials are:
  !
  !  P1(x) = 1
  !  P2(x) = x
  !  P3(x) = 1/2 x^2 - 1/6
  !
  ! The first one is obviously the constant polynomial and the second and thrid
  ! have vanishing mean-value on the reference element [-1,1].

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T2_1D-element always computes function value and 1st derivative.
  ! That is even faster than when using two IF commands for preventing
  ! the computation of one of the values!

  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC1D)) then
    Dbas(1,DER_FUNC1D) = 1.0_DP
    Dbas(2,DER_FUNC1D) = Dpoint(1)
    Dbas(3,DER_FUNC1D) = 0.5_DP*Dpoint(1)*Dpoint(1)-1.0_DP/6.0_DP
!  endif

  ! If x-derivatives are desired, calculate them.
!  if (Bder(DER_DERIV1D_X)) then
    Dbas(1,DER_DERIV1D_X) = 0.0_DP
    Dbas(2,DER_DERIV1D_X) = 1.0_DP / ddetj
    Dbas(3,DER_DERIV1D_X) = Dpoint(1) / ddetj
!  endif

  ! If xx-derivatives are desired, calculate them.
!  if (Bder(DER_DERIV1D_XX)) then
    Dbas(1,DER_DERIV1D_XX) = 0.0_DP
    Dbas(2,DER_DERIV1D_XX) = 0.0_DP
    Dbas(3,DER_DERIV1D_XX) = 1.0_DP / ddetj
!  endif


  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine elem_DG_T2_1D_mult (celement, Dcoords, Djac, Ddetj, &
                                      Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T2_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
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

  integer :: i   ! point counter

  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The DG_T2_1D-element always computes function value and 1st derivatives
  ! That is even faster than when using two IF commands for preventing
  ! the computation of one of the values!

  !if function values are desired
  !IF (Bder(DER_FUNC1D)) THEN
    do i=1,npoints
      Dbas(1,DER_FUNC1D,i) = 1.0_DP
      Dbas(2,DER_FUNC1D,i) = Dpoints(1,i)
      Dbas(3,DER_FUNC1D,i) = 0.5_DP*Dpoints(1,i)*Dpoints(1,i)-1.0_DP/6.0_DP
    !end do
  !ENDIF

  !if x-derivatives are desired
!  IF (Bder(DER_DERIV1D_X)) THEN
    !do i=1,npoints
      Dbas(1,DER_DERIV1D_X,i) = 0.0_DP
      Dbas(2,DER_DERIV1D_X,i) = 1.0_DP / Ddetj(i)
      Dbas(3,DER_DERIV1D_X,i) = Dpoints(1,i) / Ddetj(i)
    !end do
!  ENDIF


  ! If xx-derivatives are desired, calculate them.
!  if (Bder(DER_DERIV1D_XX)) then
    !do i=1,npoints
      Dbas(1,DER_DERIV1D_XX,i) = 0.0_DP
      Dbas(2,DER_DERIV1D_XX,i) = 0.0_DP
      Dbas(3,DER_DERIV1D_XX,i) = 1.0_DP / Ddetj(i)
   end do
!  endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_DG_T2_1D_sim (celement, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, nelements, &
                                Dpoints, rperfconfig)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_DG_T2_1D.
  integer(I32), intent(in)  :: celement

  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints

  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
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

  integer :: i   ! point counter
  integer :: j   ! element counter

    !if function values are desired
    if (Bder(DER_FUNC1D)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_FUNC1D,i,j) = 1.0_DP
          Dbas(2,DER_FUNC1D,i,j) = Dpoints(1,i,j)
          Dbas(3,DER_FUNC1D,i,j) = 0.5_DP*Dpoints(1,i,j)*Dpoints(1,i,j)-1.0_DP/6.0_DP
        end do

      end do
      !$omp end parallel do

    end if

    !if x-derivatives are desired
    if (Bder(DER_DERIV1D_X)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_DERIV1D_X,i,j) = 0.0_DP
          Dbas(2,DER_DERIV1D_X,i,j) = 1.0_DP / Ddetj(i,j)
          Dbas(3,DER_DERIV1D_X,i,j) = Dpoints(1,i,j) / Ddetj(i,j)
        end do

      end do
      !$omp end parallel do

    end if

    !if xx-derivatives are desired
    if (Bder(DER_DERIV1D_XX)) then

      !$omp parallel do default(shared) private(i) &
      !$omp if(nelements > rperfconfig%NELEMMIN_OMP)
      do j=1,nelements

        do i=1,npoints
          Dbas(1,DER_DERIV1D_XX,i,j) = 0.0_DP
          Dbas(2,DER_DERIV1D_XX,i,j) = 0.0_DP / Ddetj(i,j)
          Dbas(3,DER_DERIV1D_XX,i,j) = 1.0_DP / Ddetj(i,j)
        end do

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

  subroutine elem_eval_P1_1D (celement, reval, Bder, Dbas)

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
  ! The P1_1D element is specified by two polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,2:
  ! {
  !   For all j = 1,...,2:
  !   {
  !     Pi(vj) = kronecker(i,j)
  !   }
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the line
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  ! P1(x) = 1/2 * (1 - x)
  ! P2(x) = 1/2 * (1 + x)

  ! Local variables
  real(DP) :: ddet,dx
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC1D,i,j) = 0.5_DP*(1.0_DP-dx)
          Dbas(2,DER_FUNC1D,i,j) = 0.5_DP*(1.0_DP+dx)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV1D_X)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,ddet) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          Dbas(1,DER_DERIV1D_X,i,j) = -0.5_DP*ddet
          Dbas(2,DER_DERIV1D_X,i,j) =  0.5_DP*ddet

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

  subroutine elem_eval_P2_1D (celement, reval, Bder, Dbas)

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
  ! The P2_1D element is specified by three polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, x^2 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,3:
  ! {
  !   For all j = 1,...,2:
  !   {
  !     Pi(vj) = kronecker(i,j)
  !   }
  !   Pi(0) = kronecker(i,3)
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the line
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  ! P1(x) = 1/2 * x * (1 - x)
  ! P2(x) = 1/2 * x * (1 + x)
  ! P3(x) = 1 - x^2

  ! Local variables
  real(DP) :: ddet,dx
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC1D,i,j) = 0.5_DP*dx*(1.0_DP-dx)
          Dbas(2,DER_FUNC1D,i,j) = 0.5_DP*dx*(1.0_DP+dx)
          Dbas(3,DER_FUNC1D,i,j) = 1.0_DP - dx*dx

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV1D_X)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,ddet) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          Dbas(1,DER_DERIV1D_X,i,j) = (dx - 0.5_DP)*ddet
          Dbas(2,DER_DERIV1D_X,i,j) = (dx + 0.5_DP)*ddet
          Dbas(3,DER_DERIV1D_X,i,j) = -2.0_DP*dx*ddet

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

  subroutine elem_eval_S31_1D (celement, reval, Bder, Dbas)

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
  ! The S31_1D element is specified by four polynomials per element.
  !
  ! The basis polynomials are constructed from the following set of monomials:
  !
  ! { 1, x, x^2, x^3 }
  !
  ! The basis polynomials Pi are constructed such that they fulfill the
  ! following conditions:
  !
  ! For all i = 1,...,4:
  ! {
  !   For all j = 1,...,2:
  !   {
  !     Pi  (vj) = kronecker(i,j)
  !     Pi_x(vj) = kronecker(i,j+2)
  !   }
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the line
  ! Pi_x being the first derivative of Pi
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  !  P1(x) = 1/4 * x * (x^2 - 3) + 1/2
  !  P2(x) = 1/4 * x * (3 - x^2) + 1/2
  !  P3(x) = 1/4 * (x + 1) * (x - 1)^2
  !  P4(x) = 1/4 * (x - 1) * (x + 1)^2
  !

  ! Local variables
  real(DP) :: ddet,dx,dlen
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(ddet,dlen,dx,i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Calculate length of element
        dlen = 0.125_DP * (reval%p_Dcoords(1,2,j) - reval%p_Dcoords(1,1,j))

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Get the determinant for this point
          ddet = reval%p_Ddetj(i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC1D,i,j) =  0.25_DP*dx*(dx*dx - 3.0_DP) + 0.5_DP
          Dbas(2,DER_FUNC1D,i,j) = -0.25_DP*dx*(dx*dx - 3.0_DP) + 0.5_DP
          Dbas(3,DER_FUNC1D,i,j) = dlen*(dx + 1.0_DP)*(dx - 1.0_DP)**2
          Dbas(4,DER_FUNC1D,i,j) = dlen*(dx - 1.0_DP)*(dx + 1.0_DP)**2

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV1D_X)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(ddet,dlen,dx,i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Calculate length of element
        dlen = 0.5_DP * (reval%p_Dcoords(1,2,j) - reval%p_Dcoords(1,1,j))

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          Dbas(1,DER_DERIV1D_X,i,j) =  0.75_DP*(dx*dx - 1.0_DP)*ddet
          Dbas(2,DER_DERIV1D_X,i,j) = -0.75_DP*(dx*dx - 1.0_DP)*ddet
          Dbas(3,DER_DERIV1D_X,i,j) = dlen*(dx*(3.0_DP*dx - 1.0_DP) - 1.0_DP)*ddet
          Dbas(4,DER_DERIV1D_X,i,j) = dlen*(dx*(3.0_DP*dx + 1.0_DP) - 1.0_DP)*ddet

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

  subroutine elem_eval_PN_1D (celement, reval, Bder, Dbas)

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
  ! TODO

  ! Local variables
  real(DP) :: ddet,dx,dx2
  integer :: i,j,k,n

    ! Get the degree of the element
    n = iand(int(ishft(celement,-16)),255)+1

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(ddet,dx,dx2,i,k)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Evaluate P1 basis functions
          Dbas(1,DER_FUNC1D,i,j) = 0.5_DP*(1.0_DP-dx)
          Dbas(2,DER_FUNC1D,i,j) = 0.5_DP*(1.0_DP+dx)

          ! Evaluate basis functions of even degree >= 2
          dx2 = 1.0_DP
          do k = 3, n, 2
            dx2 = dx2*dx*dx
            Dbas(k,DER_FUNC1D,i,j) = 1.0_DP - dx2
          end do

          ! Evaluate basis functions of odd degree >= 3
          dx2 = 1.0_DP
          do k = 4, n, 2
            dx2 = dx2*dx*dx
            Dbas(k,DER_FUNC1D,i,j) = (dx2 - 1.0_DP)*dx
          end do

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV1D_X)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(ddet,dx,dx2,i,k)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element of P1 basis functions
          Dbas(1,DER_DERIV1D_X,i,j) = -0.5_DP*ddet
          Dbas(2,DER_DERIV1D_X,i,j) =  0.5_DP*ddet

          ! X-derivatives on real element of basis functions of
          ! even degree >= 2
          dx2 = dx*ddet
          do k = 3, n, 2
            Dbas(k,DER_DERIV1D_X,i,j) = -real(k-1,DP)*dx2
            dx2 = dx2*dx*dx
          end do

          ! X-derivatives on real element of basis functions of
          ! odd degree >= 3
          dx2 = ddet
          do k = 4, n, 2
            dx2 = dx2*dx*dx
            Dbas(k,DER_DERIV1D_X,i,j) = real(k-1,DP)*dx2 - ddet
          end do

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

end module
