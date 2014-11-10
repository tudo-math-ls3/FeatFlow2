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

!$ use omp_lib
  use fsystem
  use elementbase
  use derivatives
  use perfconfig

  implicit none

  private

  public :: elem_eval_P0_3D
  public :: elem_eval_P1_3D

contains

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

 subroutine elem_eval_P0_3D (celement, reval, Bder, Dbas)

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
  ! The P0_3D element is specified by a single basis function, constant in the element.
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

 subroutine elem_eval_P1_3D (celement, reval, Bder, Dbas)

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
  real(DP) :: ddet
  integer :: i,j

  ! A matrix for the inverse Jacobian matrix
  real(DP), dimension(9) :: Dinv

    ! Calculate function values?
    if(Bder(DER_FUNC3D)) then

      ! Loop through all elements
      !$omp parallel do default(shared) private(i)&
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j = 1, reval%nelements

        ! Loop through all points on the current element
        do i = 1, reval%npointsPerElement

          ! Evaluate basis functions
          Dbas(1,DER_FUNC3D,i,j) = reval%p_DpointsRef(1,i,j)
          Dbas(2,DER_FUNC3D,i,j) = reval%p_DpointsRef(2,i,j)
          Dbas(3,DER_FUNC3D,i,j) = reval%p_DpointsRef(3,i,j)
          Dbas(4,DER_FUNC3D,i,j) = reval%p_DpointsRef(4,i,j)
          
        end do ! i

      end do ! j
      !$omp end parallel do

    end if

    ! Calculate derivatives?
    if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
        (Bder(DER_DERIV3D_Z))) then
      
      !$omp parallel do default(shared) private(i,ddet,Dinv) &
      !$omp if(reval%nelements > reval%p_rperfconfig%NELEMMIN_OMP)
      do j=1,reval%nelements
        
        ! Since the Jacobian matrix (and therefore also its determinant) is
        ! constant for all points on one element, we only need to invert
        ! the matrix once per element.
        ddet = 1.0_DP / reval%p_Ddetj(1,j)
        Dinv(1)=(reval%p_Djac(5,1,j)*reval%p_Djac(9,1,j)-reval%p_Djac(8,1,j)*reval%p_Djac(6,1,j))*ddet
        Dinv(2)=(reval%p_Djac(8,1,j)*reval%p_Djac(3,1,j)-reval%p_Djac(2,1,j)*reval%p_Djac(9,1,j))*ddet
        Dinv(3)=(reval%p_Djac(2,1,j)*reval%p_Djac(6,1,j)-reval%p_Djac(5,1,j)*reval%p_Djac(3,1,j))*ddet
        Dinv(4)=(reval%p_Djac(7,1,j)*reval%p_Djac(6,1,j)-reval%p_Djac(4,1,j)*reval%p_Djac(9,1,j))*ddet
        Dinv(5)=(reval%p_Djac(1,1,j)*reval%p_Djac(9,1,j)-reval%p_Djac(7,1,j)*reval%p_Djac(3,1,j))*ddet
        Dinv(6)=(reval%p_Djac(4,1,j)*reval%p_Djac(3,1,j)-reval%p_Djac(1,1,j)*reval%p_Djac(6,1,j))*ddet
        Dinv(7)=(reval%p_Djac(4,1,j)*reval%p_Djac(8,1,j)-reval%p_Djac(7,1,j)*reval%p_Djac(5,1,j))*ddet
        Dinv(8)=(reval%p_Djac(7,1,j)*reval%p_Djac(2,1,j)-reval%p_Djac(1,1,j)*reval%p_Djac(8,1,j))*ddet
        Dinv(9)=(reval%p_Djac(1,1,j)*reval%p_Djac(5,1,j)-reval%p_Djac(4,1,j)*reval%p_Djac(2,1,j))*ddet
        
        ! Loop through all points on the current element
        do i=1,reval%npointsPerElement

          ! X-derivatives on real element
          Dbas(1,DER_DERIV3D_X,i,j) = -(Dinv(1)+Dinv(2)+Dinv(3))
          Dbas(2,DER_DERIV3D_X,i,j) = Dinv(1)
          Dbas(3,DER_DERIV3D_X,i,j) = Dinv(2)
          Dbas(4,DER_DERIV3D_X,i,j) = Dinv(3)
          
          ! Y-derivatives on real element
          Dbas(1,DER_DERIV3D_Y,i,j) = -(Dinv(4)+Dinv(5)+Dinv(6))
          Dbas(2,DER_DERIV3D_Y,i,j) = Dinv(4)
          Dbas(3,DER_DERIV3D_Y,i,j) = Dinv(5)
          Dbas(4,DER_DERIV3D_Y,i,j) = Dinv(6)
          
          ! Z-derivatives on real element
          Dbas(1,DER_DERIV3D_Z,i,j) = -(Dinv(7)+Dinv(8)+Dinv(9))
          Dbas(2,DER_DERIV3D_Z,i,j) = Dinv(7)
          Dbas(3,DER_DERIV3D_Z,i,j) = Dinv(8)
          Dbas(4,DER_DERIV3D_Z,i,j) = Dinv(9)
        end do

      end do
      !$omp end parallel do
      
    end if
    
  end subroutine

end module
