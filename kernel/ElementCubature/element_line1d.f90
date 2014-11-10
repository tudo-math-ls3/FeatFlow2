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

  public :: elem_eval_P0_1D
  public :: elem_eval_P1_1D
  public :: elem_eval_P2_1D
  public :: elem_eval_S31_1D
  public :: elem_eval_PN_1D
  public :: elem_eval_DG_T0_1D
  public :: elem_eval_DG_T1_1D
  public :: elem_eval_DG_T2_1D

contains

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

  subroutine elem_eval_P0_1D (celement, reval, Bder, Dbas)

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
  ! The P0_1D element is specified by a single basis function per element.
  !
  ! The function value of the basis function is =1, the derivatives are all 0!

  ! Local variables
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas(1,DER_FUNC1D,i,j) = 1.0_DP

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

          ! Get Jacobian determinant
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
  !   Pi(x) = kronecker(i,3)
  ! }
  !
  ! With:
  ! vj being the j-th local corner vertice of the line
  ! x being the midpoint of the line
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  ! P1(x) = 1/2 * x * (x - 1)
  ! P2(x) = 1/2 * x * (x + 1)
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
          Dbas(1,DER_FUNC1D,i,j) = 0.5_DP*dx*(dx-1.0_DP)
          Dbas(2,DER_FUNC1D,i,j) = 0.5_DP*dx*(dx+1.0_DP)
          Dbas(3,DER_FUNC1D,i,j) = (1.0_DP - dx)*(1.0_DP + dx)

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

          ! Get Jacobian determinant
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
  ! Pi_x being the first derivative of Pi in real coordinates
  !
  ! On the reference element, the above combination of monomial set and
  ! basis polynomial conditions leads to the following basis polynomials:
  !
  !  P1(x) = 1/4 * x * (x^2 - 3) + 1/2
  !  P2(x) = 1/4 * x * (3 - x^2) + 1/2
  !  P3(x) = 1/4 * L * (x + 1) * (x - 1)^2
  !  P4(x) = 1/4 * L * (x - 1) * (x + 1)^2
  !
  ! where L := 1/2 * (v_{i+1} - v_{i}) denotes the half length of the
  ! element.

  ! Local variables
  real(DP) :: ddet,dx,dlen
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(dlen,dx,i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Calculate length of element
        dlen = 0.5_DP * (reval%p_Dcoords(1,2,j) - reval%p_Dcoords(1,1,j))

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC1D,i,j) =  0.25_DP*dx*(dx*dx - 3.0_DP) + 0.5_DP
          Dbas(2,DER_FUNC1D,i,j) = -0.25_DP*dx*(dx*dx - 3.0_DP) + 0.5_DP
          Dbas(3,DER_FUNC1D,i,j) =  0.25_DP*dlen*(dx + 1.0_DP)*(dx - 1.0_DP)**2
          Dbas(4,DER_FUNC1D,i,j) =  0.25_DP*dlen*(dx - 1.0_DP)*(dx + 1.0_DP)**2

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

          ! Get Jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          Dbas(1,DER_DERIV1D_X,i,j) =  0.75_DP*(dx - 1.0_DP)*(dx + 1.0_DP)*ddet
          Dbas(2,DER_DERIV1D_X,i,j) = -0.75_DP*(dx - 1.0_DP)*(dx + 1.0_DP)*ddet
          Dbas(3,DER_DERIV1D_X,i,j) =  0.25_DP*dlen*(3.0_DP*dx + 1.0_DP)*(dx - 1.0_DP)*ddet
          Dbas(4,DER_DERIV1D_X,i,j) =  0.25_DP*dlen*(3.0_DP*dx - 1.0_DP)*(dx + 1.0_DP)*ddet

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

          ! Get Jacobian determinant
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

  !************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine elem_eval_DG_T0_1D (celement, reval, Bder, Dbas)

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
  ! The P0_1D element is specified by a single basis function per element.
  !
  ! The function value of the basis function is =1, the derivatives are all 0!

  ! Local variables
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas(1,DER_FUNC1D,i,j) = 1.0_DP

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

  subroutine elem_eval_DG_T1_1D (celement, reval, Bder, Dbas)

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
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC1D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas(1,DER_FUNC1D,i,j) = 1.0_DP
          Dbas(2,DER_FUNC1D,i,j) = reval%p_DpointsRef(1,i,j)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV1D_X)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! X-derivatives on real element
          Dbas(1,DER_DERIV1D_X,i,j) = 0.0_DP
          Dbas(2,DER_DERIV1D_X,i,j) = 1.0_DP / reval%p_Ddetj(i,j)

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

  subroutine elem_eval_DG_T2_1D (celement, reval, Bder, Dbas)

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
          Dbas(1,DER_FUNC1D,i,j) = 1.0_DP
          Dbas(2,DER_FUNC1D,i,j) = dx
          Dbas(3,DER_FUNC1D,i,j) = 0.5_DP*dx*dx-1.0_DP/6.0_DP

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

          ! Get Jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          Dbas(1,DER_DERIV1D_X,i,j) = 0.0_DP
          Dbas(2,DER_DERIV1D_X,i,j) = ddet
          Dbas(3,DER_DERIV1D_X,i,j) = dx*ddet

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV1D_XX)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,ddet) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get Jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          Dbas(1,DER_DERIV1D_X,i,j) = 0.0_DP
          Dbas(2,DER_DERIV1D_X,i,j) = 0.0_DP
          Dbas(3,DER_DERIV1D_X,i,j) = ddet

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

end module
