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

!$ use omp_lib
  use basicgeometry
  use derivatives
  use elementbase
  use fsystem
  use perfconfig

  implicit none

  private

  public :: elem_eval_R0_3D
  public :: elem_eval_R1_3D

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

  subroutine elem_eval_R0_3D (celement, reval, Bder, Dbas)

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
  ! The R0_3D element is specified by a single basis function, constant in the element.
  !
  ! The function value of the basis function is =1, the derivatives are all 0!

  ! Local variables
  integer :: i,j

    ! Calculate function values?
    if(Bder(DER_FUNC3D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas(1,DER_FUNC3D,i,j) = 1.0_DP

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

  subroutine elem_eval_R1_3D (celement, reval, Bder, Dbas)

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
  ! The R1_3D element is specified by six polynomials on the reference element.
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

  ! Local variables
  real(DP) ::  djx,djy,djz,ddet
  integer :: i,j

  !auxiliary vector containing the first derivatives on the reference element
  real(DP), dimension(6,NDIM3D) :: Dhelp

    ! Calculate function values?
    if(Bder(DER_FUNC3D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas(1,DER_FUNC3D,i,j) = 0.25_DP*(1.0_DP-reval%p_DpointsRef(1,i,j)-reval%p_DpointsRef(2,i,j))&
                                          *(1.0_DP-reval%p_DpointsRef(3,i,j))
          Dbas(2,DER_FUNC3D,i,j) = 0.25_DP*reval%p_DpointsRef(1,i,j)*(1.0_DP-reval%p_DpointsRef(3,i,j))
          Dbas(3,DER_FUNC3D,i,j) = 0.25_DP*reval%p_DpointsRef(2,i,j)*(1.0_DP-reval%p_DpointsRef(3,i,j))
          Dbas(4,DER_FUNC3D,i,j) = 0.25_DP*(1.0_DP-reval%p_DpointsRef(1,i,j)-reval%p_DpointsRef(2,i,j))&
                                          *(1.0_DP+reval%p_DpointsRef(3,i,j))
          Dbas(5,DER_FUNC3D,i,j) = 0.25_DP*reval%p_DpointsRef(1,i,j)*(1.0_DP+reval%p_DpointsRef(3,i,j))
          Dbas(6,DER_FUNC3D,i,j) = 0.25_DP*reval%p_DpointsRef(2,i,j)*(1.0_DP+reval%p_DpointsRef(3,i,j))

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
        (Bder(DER_DERIV3D_Z))) then
      
      !$omp parallel do default(shared) private(i,ddet,djx,djy,djz,Dhelp) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j=1,reval%nelements

        ! Loop through all points on the current element
        do i=1,reval%npointsPerElement

          Dhelp(1,1) =-(1.0_DP-reval%p_DpointsRef(3,i,j))
          Dhelp(2,1) = (1.0_DP-reval%p_DpointsRef(3,i,j))
          Dhelp(3,1) = 0.0_DP
          Dhelp(4,1) =-(1.0_DP+reval%p_DpointsRef(3,i,j))
          Dhelp(5,1) = (1.0_DP+reval%p_DpointsRef(3,i,j))
          Dhelp(6,1) = 0.0_DP
          Dhelp(1,2) =-(1.0_DP-reval%p_DpointsRef(3,i,j))
          Dhelp(2,2) = 0.0_DP
          Dhelp(3,2) = (1.0_DP-reval%p_DpointsRef(3,i,j))
          Dhelp(4,2) =-(1.0_DP+reval%p_DpointsRef(3,i,j))
          Dhelp(5,2) = 0.0_DP
          Dhelp(6,2) = (1.0_DP+reval%p_DpointsRef(3,i,j))
          Dhelp(1,3) =-(1.0_DP-reval%p_DpointsRef(1,i,j)-reval%p_DpointsRef(2,i,j))
          Dhelp(2,3) =-reval%p_DpointsRef(1,i,j)
          Dhelp(3,3) =-reval%p_DpointsRef(2,i,j)
          Dhelp(4,3) =+(1.0_DP-reval%p_DpointsRef(1,i,j)-reval%p_DpointsRef(2,i,j))
          Dhelp(5,3) =+reval%p_DpointsRef(1,i,j)
          Dhelp(6,3) =+reval%p_DpointsRef(2,i,j)
          
          ddet = 0.25_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          djx = reval%p_Djac(5,i,j)*reval%p_Djac(9,i,j) - reval%p_Djac(6,i,j)*reval%p_Djac(8,i,j)
          djy = reval%p_Djac(8,i,j)*reval%p_Djac(3,i,j) - reval%p_Djac(2,i,j)*reval%p_Djac(9,i,j)
          djz = reval%p_Djac(2,i,j)*reval%p_Djac(6,i,j) - reval%p_Djac(5,i,j)*reval%p_Djac(3,i,j)
          Dbas(1,DER_DERIV3D_X,i,j) = ddet * &
              (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
          Dbas(2,DER_DERIV3D_X,i,j) = ddet * &
              (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
          Dbas(3,DER_DERIV3D_X,i,j) = ddet * &
              (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
          Dbas(4,DER_DERIV3D_X,i,j) = ddet * &
              (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
          Dbas(5,DER_DERIV3D_X,i,j) = ddet * &
              (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
          Dbas(6,DER_DERIV3D_X,i,j) = ddet * &
              (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
          
          ! Y-derivatives on real element
          djx = reval%p_Djac(7,i,j)*reval%p_Djac(6,i,j) - reval%p_Djac(4,i,j)*reval%p_Djac(9,i,j)
          djy = reval%p_Djac(1,i,j)*reval%p_Djac(9,i,j) - reval%p_Djac(7,i,j)*reval%p_Djac(3,i,j)
          djz = reval%p_Djac(4,i,j)*reval%p_Djac(3,i,j) - reval%p_Djac(1,i,j)*reval%p_Djac(6,i,j)
          Dbas(1,DER_DERIV3D_Y,i,j) = ddet * &
              (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
          Dbas(2,DER_DERIV3D_Y,i,j) = ddet * &
              (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
          Dbas(3,DER_DERIV3D_Y,i,j) = ddet * &
              (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
          Dbas(4,DER_DERIV3D_Y,i,j) = ddet * &
              (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
          Dbas(5,DER_DERIV3D_Y,i,j) = ddet * &
              (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
          Dbas(6,DER_DERIV3D_Y,i,j) = ddet * &
              (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
          
          ! Z-derivatives on real element
          djx = reval%p_Djac(4,i,j)*reval%p_Djac(8,i,j) - reval%p_Djac(7,i,j)*reval%p_Djac(5,i,j)
          djy = reval%p_Djac(7,i,j)*reval%p_Djac(2,i,j) - reval%p_Djac(1,i,j)*reval%p_Djac(8,i,j)
          djz = reval%p_Djac(1,i,j)*reval%p_Djac(5,i,j) - reval%p_Djac(4,i,j)*reval%p_Djac(2,i,j)
          Dbas(1,DER_DERIV3D_Z,i,j) = ddet * &
              (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
          Dbas(2,DER_DERIV3D_Z,i,j) = ddet * &
              (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
          Dbas(3,DER_DERIV3D_Z,i,j) = ddet * &
              (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
          Dbas(4,DER_DERIV3D_Z,i,j) = ddet * &
              (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
          Dbas(5,DER_DERIV3D_Z,i,j) = ddet * &
              (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
          Dbas(6,DER_DERIV3D_Z,i,j) = ddet * &
              (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
          
        end do

      end do
      !$omp end parallel do
      
    end if

  end subroutine

end module
