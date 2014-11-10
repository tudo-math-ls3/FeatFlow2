!##############################################################################
!# ****************************************************************************
!# <name> element_pyra3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the implementations of the 3D pyramid basis
!# functions.
!#
!# </purpose>
!##############################################################################

module element_pyra3d

!$ use omp_lib
  use fsystem
  use basicgeometry
  use elementbase
  use derivatives
  use transformation
  use perfconfig

  implicit none

  private

  public :: elem_eval_Y0_3D
  public :: elem_eval_Y1_3D

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

 subroutine elem_eval_Y0_3D (celement, reval, Bder, Dbas)

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
  ! The Y0_3D element is specified by a single basis function, constant in the element.
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

  subroutine elem_eval_Y1_3D (celement, reval, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
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
  ! The Y1_3D element is specified by five polynomials on the reference element.
  ! These five polynomials are:
  !
  !  P1(x,y,z) = 1/4 * (1-x) * (1-y) * (1-z)
  !  P2(x,y,z) = 1/4 * (1+x) * (1-y) * (1-z)
  !  P3(x,y,z) = 1/4 * (1+x) * (1+y) * (1-z)
  !  P4(x,y,z) = 1/4 * (1-x) * (1+y) * (1-z)
  !  P5(x,y,z) = z
  !
  ! Each of them calculated that way that Pi(Vj)=delta_ij (Kronecker)
  ! for V1,...,V5 being the five corner vertices of the reference pyramid.

  ! Parameter: number of local basis functions
  integer, parameter :: NBAS = 5

  ! Local variables
  real(DP) :: ddet,dx,dy,dz
  integer :: i,j

  ! derivatives on reference element
  real(DP), dimension(NBAS,NDIM3D) :: DrefDer

    ! Calculate function values?
    if(Bder(DER_FUNC3D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy,dz) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)
          dz = reval%p_DpointsRef(3,i,j)

          ! Evaluate basis functions
          Dbas(1,DER_FUNC3D,i,j) = 0.25_DP*(1_DP-dx)*(1_DP-dy)*(1_DP-dz)
          Dbas(2,DER_FUNC3D,i,j) = 0.25_DP*(1_DP+dx)*(1_DP-dy)*(1_DP-dz)
          Dbas(3,DER_FUNC3D,i,j) = 0.25_DP*(1_DP+dx)*(1_DP+dy)*(1_DP-dz)
          Dbas(4,DER_FUNC3D,i,j) = 0.25_DP*(1_DP-dx)*(1_DP+dy)*(1_DP-dz)
          Dbas(5,DER_FUNC3D,i,j) = dz

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if(Bder(DER_DERIV3D_X) .or. Bder(DER_DERIV3D_Y) .or. Bder(DER_DERIV3D_Z)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i,dx,dy,dz,ddet,DrefDer) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Get the point coordinates
          dx = reval%p_DpointsRef(1,i,j)
          dy = reval%p_DpointsRef(2,i,j)
          dz = reval%p_DpointsRef(3,i,j)

          ! Calculate derivatives on reference element
          ! X-derivatives
          DrefDer(1,1) =-0.25_DP*(1_DP-dy)*(1_DP-dz)
          DrefDer(2,1) = 0.25_DP*(1_DP-dy)*(1_DP-dz)
          DrefDer(3,1) = 0.25_DP*(1_DP+dy)*(1_DP-dz)
          DrefDer(4,1) =-0.25_DP*(1_DP+dy)*(1_DP-dz)
          DrefDer(5,1) = 0_DP
          ! Y-derivatives
          DrefDer(1,2) =-0.25_DP*(1_DP-dx)*(1_DP-dz)
          DrefDer(2,2) =-0.25_DP*(1_DP+dx)*(1_DP-dz)
          DrefDer(3,2) = 0.25_DP*(1_DP+dx)*(1_DP-dz)
          DrefDer(4,2) = 0.25_DP*(1_DP-dx)*(1_DP-dz)
          DrefDer(5,2) = 0_DP
          ! Z-derivatives
          DrefDer(1,3) =-0.25_DP*(1_DP-dx)*(1_DP-dy)
          DrefDer(2,3) =-0.25_DP*(1_DP+dx)*(1_DP-dy)
          DrefDer(3,3) =-0.25_DP*(1_DP+dx)*(1_DP+dy)
          DrefDer(4,3) =-0.25_DP*(1_DP-dx)*(1_DP+dy)
          DrefDer(5,3) = 1_DP

          ! Remark: Please note that the following code is universal and does
          ! not need to be modified for other parametric 3D pyramid elements!

          ! Get jacobian determinant
          ddet = 1.0_DP / reval%p_Ddetj(i,j)

          ! X-derivatives on real element
          dx = (reval%p_Djac(5,i,j)*reval%p_Djac(9,i,j)&
               -reval%p_Djac(6,i,j)*reval%p_Djac(8,i,j))*ddet
          dy = (reval%p_Djac(8,i,j)*reval%p_Djac(3,i,j)&
               -reval%p_Djac(2,i,j)*reval%p_Djac(9,i,j))*ddet
          dz = (reval%p_Djac(2,i,j)*reval%p_Djac(6,i,j)&
               -reval%p_Djac(5,i,j)*reval%p_Djac(3,i,j))*ddet
          Dbas(1:NBAS,DER_DERIV3D_X,i,j) = dx*DrefDer(1:NBAS,1) &
                  + dy*DrefDer(1:NBAS,2) + dz*DrefDer(1:NBAS,3)

          ! Y-derivatives on real element
          dx = (reval%p_Djac(7,i,j)*reval%p_Djac(6,i,j)&
               -reval%p_Djac(4,i,j)*reval%p_Djac(9,i,j))*ddet
          dy = (reval%p_Djac(1,i,j)*reval%p_Djac(9,i,j)&
               -reval%p_Djac(7,i,j)*reval%p_Djac(3,i,j))*ddet
          dz = (reval%p_Djac(4,i,j)*reval%p_Djac(3,i,j)&
               -reval%p_Djac(1,i,j)*reval%p_Djac(6,i,j))*ddet
          Dbas(1:NBAS,DER_DERIV3D_Y,i,j) = dx*DrefDer(1:NBAS,1) &
                  + dy*DrefDer(1:NBAS,2) + dz*DrefDer(1:NBAS,3)

          ! Z-derivatives on real element
          dx = (reval%p_Djac(4,i,j)*reval%p_Djac(8,i,j)&
               -reval%p_Djac(7,i,j)*reval%p_Djac(5,i,j))*ddet
          dy = (reval%p_Djac(7,i,j)*reval%p_Djac(2,i,j)&
               -reval%p_Djac(1,i,j)*reval%p_Djac(8,i,j))*ddet
          dz = (reval%p_Djac(1,i,j)*reval%p_Djac(5,i,j)&
               -reval%p_Djac(4,i,j)*reval%p_Djac(2,i,j))*ddet
          Dbas(1:NBAS,DER_DERIV3D_Z,i,j) = dx*DrefDer(1:NBAS,1) &
                  + dy*DrefDer(1:NBAS,2) + dz*DrefDer(1:NBAS,3)

        end do ! i

      end do ! j
      !$omp end parallel do

    end if

  end subroutine

end module
