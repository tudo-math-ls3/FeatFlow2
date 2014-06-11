module ext_visco

use fsystem
use storage
use cubature
use derivatives
use triangulation
use element
use spatialdiscretisation
use linearsystemscalar
use linearsystemblock
use blockmatassemblybase
use blockmatassembly
use feevaluation2
use collection, only: t_collection

implicit none

contains

  ! Function: nu_poliquen(u_n, p_n)
  real(DP) function func_nu_poliquen(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps) result(nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps)
    nu = sqrt(0.5_DP) * nu0 * ((p + deps)/y + (p + deps)/(sqrt((p + deps)) + y))

  end function

  ! Function: d_1 nu_poliquen(u_n, p_n)
  real(DP) function func_d1nu_poliquen(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps) result(d1nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps)
    d1nu = -sqrt(0.125_DP) * nu0 * (p + deps)/y * ( 1/(y*y) + 1/(sqrt((p + deps)) + y)**2 )

  end function

  ! Function: d_2 nu_poliquen(u_n, p_n)
  real(DP) function func_d2nu_poliquen(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps) result(d2nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps)
    d2nu = sqrt(0.5_DP) * nu0 * ( 1/y + (0.5_DP*sqrt((p + deps)) + y)/(sqrt((p + deps))+y)**2 )

  end function

  !****************************************************************************

  ! Function: nu_power(u_n, p_n)
  real(DP) function func_nu_power(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps) result(nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps) 
    nu = nu0 * (y**(dpower - 2.0_DP))

  end function

  ! Function: d_1 nu_power(u_n, p_n)
  real(DP) function func_d1nu_power(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps) result(d1nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps)
    d1nu = (dpower/2.0_DP - 1.0_DP) * nu0 * y**(dpower - 4)

  end function

  ! Function: d_2 nu_power(u_n, p_n)
  real(DP) function func_d2nu_power(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps) result(d2nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps)
    d2nu = 0.0_DP

  end function

  !****************************************************************************

  ! Function: nu_bingham(u_n, p_n)
  real(DP) function func_nu_bingham(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps) result(nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps**2) 
    nu = nu0 + 0.5_DP * sqrt(2.0) * dyield * y**(-1.0_DP)

  end function

  ! Function: d_1 nu_bingham(u_n, p_n)
  real(DP) function func_d1nu_bingham(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps) result(d1nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps**2)
    d1nu = -0.25_DP * sqrt(2.0) * dyield * y**(-3.0_DP)

  end function

  ! Function: d_2 nu_bingham(u_n, p_n)
  real(DP) function func_d2nu_bingham(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps) result(d2nu)
  real(DP), intent(in) :: u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps
  real(DP) :: y

    y = sqrt(dxu1**2 + 0.5_DP * (dyu1 + dxu2)**2 + dyu2**2 + deps**2)
    d2nu = 0.0_DP

  end function

  !****************************************************************************



!<subroutine>

  subroutine fcalc_extVisco(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

  ! Local variables
  real(DP) :: nu0, nu, d1nu, d2nu, dxpsi, dypsi, dxphi, dyphi, phi, psi
  real(DP) :: u1, u2, p, dxu1, dyu1, dxu2, dyu2
  real(DP) :: dsl0, dsl1, dsl2, dsk0, dsk1, dsk2, dsn0, dsn1, dsm0, dsb0, dsd0
  real(DP) :: dpower, deps, dyield
  integer :: iel, icubp, idofe, jdofe,ndofv,ndofp
  real(DP), dimension(:,:), pointer :: p_Domega
  real(DP), dimension(:,:,:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,p_DD1,p_DD2, p_Dp
  real(DP), dimension(:,:,:,:), pointer :: p_Dvelo, p_Dpres, p_Du

    p_Domega => rassemblyData%p_DcubWeight
    p_DA11 => RmatrixData(1,1)%p_Dentry
    p_DA12 => RmatrixData(1,2)%p_Dentry
    p_DA21 => RmatrixData(2,1)%p_Dentry
    p_DA22 => RmatrixData(2,2)%p_Dentry
    p_DB1 => RmatrixData(1,3)%p_Dentry
    p_DB2 => RmatrixData(2,3)%p_Dentry
    p_DD1 => RmatrixData(3,1)%p_Dentry
    p_DD2 => RmatrixData(3,2)%p_Dentry
    p_Dvelo => RmatrixData(1,1)%p_DbasTrial
    p_Dpres => RmatrixData(1,3)%p_DbasTrial
    ndofv = RmatrixData(1,1)%ndofTrial
    ndofp = RmatrixData(1,3)%ndofTrial

    p_Du => revalVectors%p_RvectorData(1)%p_DdataVec
    p_Dp => revalVectors%p_RvectorData(2)%p_Ddata

    ! fetch scaling factors
    nu0  = rcollection%DquickAccess(1) ! nu_0
    if (rcollection%IquickAccess(2) .eq. 0) then
      ! gradient tensor
      dsl0 = rcollection%DquickAccess(2) ! scale for L
      dsl1 = rcollection%DquickAccess(3) ! scale for L*
      dsl2 = rcollection%DquickAccess(8) ! scale for B*
      dsk0 = 0.0_DP
      dsk1 = 0.0_DP
      dsk2 = 0.0_DP
    else
      ! deformation tensor
      dsk0 = rcollection%DquickAccess(2) ! scale for L
      dsk1 = rcollection%DquickAccess(3) ! scale for L*
      dsk2 = rcollection%DquickAccess(8) ! scale for B*
      dsl0 = 0.0_DP
      dsl1 = 0.0_DP
      dsl2 = 0.0_DP
    endif
    dsn0 = rcollection%DquickAccess(4) ! scale for N
    dsn1 = rcollection%DquickAccess(5) ! scale for N*
    dsm0 = rcollection%DquickAccess(6) ! scale for M
    dsb0 = rcollection%DquickAccess(7) ! scale for B
    dsd0 = rcollection%DquickAccess(9) ! scale for D

    dpower = rcollection%DquickAccess(10)    ! index for power law
    deps = rcollection%DquickAccess(11)      ! regularization parameter
    dyield = rcollection%DquickAccess(12)    ! yield for Bingham fluid
 

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! evaluate velocity and pressure
        u1 = p_Du(1,icubp,iel,DER_FUNC2D)
        u2 = p_Du(2,icubp,iel,DER_FUNC2D)
        dxu1 = p_Du(1,icubp,iel,DER_DERIV2D_X)
        dyu1 = p_Du(1,icubp,iel,DER_DERIV2D_Y)
        dxu2 = p_Du(2,icubp,iel,DER_DERIV2D_X)
        dyu2 = p_Du(2,icubp,iel,DER_DERIV2D_Y)
        p  = p_Dp(icubp,iel,DER_FUNC2D)

        ! compute nu and its derivatives
        if (rcollection%IquickAccess(1) .eq. 100) then
          nu = nu0
          d1nu = 0.0_DP
          d2nu = 0.0_DP

        else if (rcollection%IquickAccess(1) .eq. 101) then
          nu = func_nu_power(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps)
          d1nu = func_d1nu_power(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps)
          d2nu = func_d2nu_power(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dpower,deps)

        else if (rcollection%IquickAccess(1) .eq. 102) then
          nu = func_nu_bingham(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps)
          d1nu = func_d1nu_bingham(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps)
          d2nu = func_d2nu_bingham(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,dyield,deps)

        else if (rcollection%IquickAccess(1) .eq. 104) then
          nu = func_nu_poliquen(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps)
          d1nu = func_d1nu_poliquen(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps)
          d2nu = func_d2nu_poliquen(u1,dxu1,dyu1,u2,dxu2,dyu2,p,nu0,deps)

        else
          nu = 0.0_DP
          d1nu = 0.0_DP
          d2nu = 0.0_DP

        end if

        ! Loop over all velocity test functions
        do idofe = 1, RmatrixData(1,1)%ndofTest

          psi = p_Dvelo(idofe,DER_FUNC2D,icubp,iel)
          dxpsi = p_Dvelo(idofe,DER_DERIV2D_X,icubp,iel)
          dypsi = p_Dvelo(idofe,DER_DERIV2D_Y,icubp,iel)

          ! Assemble velocity blocks A
          do jdofe = 1, RmatrixData(1,1)%ndofTrial

            phi = p_Dvelo(jdofe,DER_FUNC2D,icubp,iel)
            dxphi = p_Dvelo(jdofe,DER_DERIV2D_X,icubp,iel)
            dyphi = p_Dvelo(jdofe,DER_DERIV2D_Y,icubp,iel)

            p_DA11(jdofe,idofe,iel) = p_DA11(jdofe,idofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== (GRAD-TENS) DIFFUSION =====
              + dsl0 * nu * (dxphi*dxpsi + dyphi*dypsi) &
              ! ===== (GRAD-TENS) DIFFUSION-FRECHET =====
              + dsl1 * d1nu * (dxu1*dxphi + dyu1*dyphi) * (dxu1*dxpsi + dxu1*dxpsi) &
              ! ===== (DEFO-TENS) DIFFUSION =====
              + dsk0 * nu * 1.0_DP * (2.0_DP*dxphi*dxpsi + dyphi*dypsi) &
              ! ===== (DEFO-TENS) DIFFUSION-FRECHET =====
              + dsk1 * d1nu * 2.0_DP * (dxu1*dxphi + 0.5_DP*(dyu1+dxu2)*dyphi) &
                                     * (dxu1*dxpsi + 0.5_DP*(dyu1+dxu2)*dypsi) &
              ! ===== CONVECTION =====
              + dsn0 * (u1*dxphi*psi + u2*dyphi*psi) &
              ! ===== CONVECTION-FRECHET =====
              + dsn1 * (dxu1 * phi * psi) &
              ! ===== REACTION =====
              + dsm0 * (phi * psi) &
              )

            p_DA12(jdofe,idofe,iel) = p_DA12(jdofe,idofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== (GRAD-TENS) DIFFUSION-FRECHET =====
              + dsl1 * d1nu * (dxu2*dxphi + dyu2*dyphi) * (dxu1*dxpsi + dyu1*dypsi) &
              ! ===== (DEFO-TENS) DIFFUSION =====
              + dsk0 * nu * 1.0_DP * (dxphi*dypsi) &
              ! ===== (DEFO-TENS) DIFFUSION-FRECHET =====
              + dsk1 * d1nu * 2.0_DP * (dyu2*dyphi + 0.5_DP*(dyu1+dxu2)*dxphi) &
                                     * (dxu1*dxpsi + 0.5_DP*(dyu1+dxu2)*dypsi) &
              ! ===== CONVECTION-FRECHET =====
              + dsn1 * (dyu1 * phi * psi) &
              )

            p_DA21(jdofe,idofe,iel) = p_DA21(jdofe,idofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== (GRAD-TENS) DIFFUSION-FRECHET =====
              + dsl1 * d1nu * (dxu1*dxphi + dyu1*dyphi) * (dxu2*dxphi + dyu2*dyphi) &
              ! ===== (DEFO-TENS) DIFFUSION =====
              + dsk0 * nu * 1.0_DP * (dyphi*dxpsi) &
              ! ===== (DEFO-TENS) DIFFUSION-FRECHET =====
              + dsk1 * d1nu * 2.0_DP * (dxu1*dxphi + 0.5_DP*(dyu1+dxu2)*dyphi) &
                                     * (dyu2*dypsi + 0.5_DP*(dyu1+dxu2)*dxpsi) &
              ! ===== CONVECTION-FRECHET =====
              + dsn1 * (dxu2 * phi * psi) &
              )

            p_DA22(jdofe,idofe,iel) = p_DA22(jdofe,idofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== (GRAD-TENS) DIFFUSION =====
              + dsl0 * nu * (dxphi*dxpsi + dyphi*dypsi) &
              ! ===== (GRAD-TENS) DIFFUSION-FRECHET =====
              + dsl1 * d1nu * (dxu2*dxphi + dyu2*dyphi) * (dxu2*dxpsi + dyu2*dyphi) &
              ! ===== (DEFO-TENS) DIFFUSION =====
              + dsk0 * nu * 1.0_DP * (2.0_DP*dyphi*dypsi + dxphi*dxpsi) &
              ! ===== (DEFO-TENS) DIFFUSION-FRECHET =====
              + dsk1 * d1nu * 2.0_DP * (dyu2*dyphi + 0.5_DP*(dyu1+dxu2)*dxphi) &
                                     * (dyu2*dypsi + 0.5_DP*(dyu1+dxu2)*dxpsi) &
              ! ===== CONVECTION =====
              + dsn0 * (u1*dxphi*psi + u2*dyphi*psi) &
              ! ===== CONVECTION-FRECHET =====
              + dsn1 * (dyu2 * phi * psi) &
              ! ===== REACTION =====
              + dsm0 * (phi * psi) &
              )

          end do ! jdofe

          ! Assemble gradient blocks B and divergence blocks D
          do jdofe=1,RmatrixData(1,3)%ndofTrial

            phi = p_Dpres(jdofe,DER_FUNC2D,icubp,iel)

            p_DB1(jdofe,idofe,iel) = p_DB1(jdofe,idofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== PRESSURE-GRADIENT =====
                - dsb0 * phi * dxpsi &
              ! ===== (GRAD-TENS) DIFFUSION-FRECHET =====
                + dsl2 * d2nu * (dxu1 * dxpsi + dyu1*dypsi) * dp &
              ! ===== (DEFO-TENS) DIFFUSION-FRECHET =====
                + dsk2 * d2nu * 2.0_DP * (dxu1 * dxpsi + 0.5_DP*(dyu1+dxu2)*dypsi) * dp &
              )

            p_DB2(jdofe,idofe,iel) = p_DB2(jdofe,idofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== PRESSURE-GRADIENT =====
                - dsb0 * phi * dypsi &
              ! ===== (GRAD-TENS) DIFFUSION-FRECHET =====
                + dsl2 * d2nu * (dxu2 * dxpsi + dyu2*dypsi) * dp &
              ! ===== (DEFO-TENS) DIFFUSION-FRECHET =====
                + dsk2 * d2nu * 2.0_DP * (dyu2 * dypsi + 0.5_DP*(dyu1+dxu2)*dxpsi) * dp &
              )

            p_DD1(idofe,jdofe,iel) = p_DD1(idofe,jdofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== VELOCITY-DIVERGENCE =====
                - dsd0 * phi * dxpsi &
              )

            p_DD2(idofe,jdofe,iel) = p_DD2(idofe,jdofe,iel) + p_Domega(icubp,iel) * &
              ( &
              ! ===== VELOCITY-DIVERGENCE =====
                - dsd0 * phi*dypsi &
              )

          end do ! jdofe

        end do ! idofe

      end do ! icubp

    end do ! iel

  end subroutine

end module
