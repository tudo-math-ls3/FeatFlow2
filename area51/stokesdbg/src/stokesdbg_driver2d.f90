module stokesdbg_driver2d

use stokesdbg_core

implicit none

contains

  ! ***********************************************************************************************

  subroutine stdrv_initBoundaryConditions(rproblem)
  type(t_problem), intent(inout) :: rproblem
  
  integer :: ivelo
  
    call parlst_getvalue_int(rproblem%p_rparam, '', 'VORTEX_VELO', ivelo, 0)

    select case(rproblem%idriver)
    case(2001)
      call stdbg_initQuadPoiseulleBCs(rproblem)
      
    case (2002)
      select case(ivelo)
      case (0,2,3)
        call stdbg_initNoSlipBCs(rproblem)
      case (1)
        call stdbg_initSlipBCs(rproblem)
      end select

    end select
  
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdrv_initFilterChain(rsystem)
  type(t_system), intent(inout) :: rsystem
  
    select case(rsystem%p_rproblem%idriver)
    case (2001)
      allocate(rsystem%p_RfilterChain(1))
      rsystem%p_RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    case (2002)
      allocate(rsystem%p_RfilterChain(2))
      rsystem%p_RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
      rsystem%p_RfilterChain(2)%ifilterType = FILTER_TOL20
      rsystem%p_RfilterChain(2)%itoL20component = 3
      
    end select
    
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdrv_postProcSol(rsystem)
  type(t_system), intent(inout) :: rsystem
  
    select case(rsystem%p_rproblem%idriver)
    case (2002)
      ! filter pressure mean
      call stdbg_filterPressureMean(rsystem)
    end select
  end subroutine

  ! ***********************************************************************************************

  subroutine stdrv_funcRhsVelo (rdiscretisation,rform,nelements,npointsPerElement,Dpoints, &
                         IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_linearForm), intent(in)                              :: rform
  integer, intent(in)                                         :: nelements
  integer, intent(in)                                         :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset
  type(t_collection), intent(inout), optional      :: rcollection
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

  integer :: idriver, icomp, ivelo
  real(DP) :: dnu,dalpha,dbeta,dgamma

    ! fetch coefficients from collection
    dnu = rcollection%DquickAccess(1)
    dalpha = rcollection%DquickAccess(2)
    dbeta = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    idriver = rcollection%IquickAccess(1)
    icomp = rcollection%IquickAccess(5)
    ivelo = rcollection%IquickAccess(6)

    Dcoefficients(1,:,:) = 0.0_DP

    select case(idriver)
    case (2001)
      ! Poiseulle-Flow
      select case(icomp)
      case (1)
        Dcoefficients(1,:,:) = 2.0_DP*dnu
      case (2)
        Dcoefficients(1,:,:) = 0.0_DP
      end select
      
    case (2002)
      ! add velocity laplacian
      select case(ivelo)
      case (0)
        ! no-flow
        Dcoefficients(1,:,:) = 0.0_DP
        
      case (1)
        select case(icomp)
        case (1)
          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) &
            +2.0_DP*dnu*SYS_PI**2 * sin(SYS_PI*Dpoints(1,:,:))*cos(SYS_PI*Dpoints(2,:,:))
        case (2)
          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) &
            -2.0_DP*dnu*SYS_PI**2 * sin(SYS_PI*Dpoints(2,:,:))*cos(SYS_PI*Dpoints(1,:,:))
        end select
        
      case (2)
        select case(icomp)
        case (1)
          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) &
            -4.0_DP*dnu*SYS_PI**2*sin(2.0_DP*SYS_PI*Dpoints(2,:,:))*&
            (2.0_DP*cos(2.0_DP*SYS_PI*Dpoints(1,:,:)) - 1.0_DP)
        case (2)
          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) &
            +4.0_DP*dnu*SYS_PI**2*sin(2.0_DP*SYS_PI*Dpoints(1,:,:))*&
            (2.0_DP*cos(2.0_DP*SYS_PI*Dpoints(2,:,:)) - 1.0_DP)
        end select
        
      case (3)
        select case(icomp)
        case (1)
          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) &
            +(2.0_DP-12.0_DP*Dpoints(1,:,:)+12.0_DP*Dpoints(1,:,:)**2)*&
            (2.0_DP*Dpoints(2,:,:)-6.0_DP*Dpoints(2,:,:)**2+4.0_DP*Dpoints(2,:,:)**3)&
            + Dpoints(1,:,:)**2*(1.0_DP-Dpoints(1,:,:))**2*(-12.0_DP+24.0_DP*Dpoints(2,:,:))
        case (2)
          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) &
            -(2.0_DP-12.0_DP*Dpoints(2,:,:)+12.0_DP*Dpoints(2,:,:)**2)*&
            (2.0_DP*Dpoints(1,:,:)-6.0_DP*Dpoints(1,:,:)**2+4.0_DP*Dpoints(1,:,:)**3)&
            - Dpoints(2,:,:)**2*(1.0_DP-Dpoints(2,:,:))**2*(-12.0_DP+24.0_DP*Dpoints(1,:,:))
        end select
      end select
    end select
    
  end subroutine

  ! ***********************************************************************************************

  subroutine stdrv_funcRhsPres (rdiscretisation,rform,nelements,npointsPerElement,Dpoints, &
                         IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_linearForm), intent(in)                              :: rform
  integer, intent(in)                                         :: nelements
  integer, intent(in)                                         :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset
  type(t_collection), intent(inout), optional      :: rcollection
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

  integer :: idriver, icomp, ipres
  real(DP) :: dnu,dalpha,dbeta,dgamma

    ! fetch coefficients from collection
    dnu = rcollection%DquickAccess(1)
    dalpha = rcollection%DquickAccess(2)
    dbeta = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    idriver = rcollection%IquickAccess(1)
    icomp = rcollection%IquickAccess(5)
    ipres = rcollection%IquickAccess(7)

    select case(idriver)
    case (2001)
      ! Poiseulle-Flow
      select case(icomp)
      case (1)
        Dcoefficients(1,:,:) = -2.0_DP*dnu*dalpha
      case (2)
        Dcoefficients(1,:,:) = 0.0_DP
      end select
      
    case (2002)
      ! add pressure gradient
      select case(ipres)
      case (0)
        Dcoefficients(1,:,:) = 0.0_DP

      case (1)
        select case(icomp)
        case (1)
          Dcoefficients(1,:,:) = -dalpha*SYS_PI*sin(SYS_PI*Dpoints(1,:,:))*cos(SYS_PI*Dpoints(2,:,:))
        case (2)
          Dcoefficients(1,:,:) = -dalpha*SYS_PI*cos(SYS_PI*Dpoints(1,:,:))*sin(SYS_PI*Dpoints(2,:,:))
        end select
        
      case (2)
        select case(icomp)
        case (1)
          Dcoefficients(1,:,:) = -dalpha*SYS_PI*cos(SYS_PI*Dpoints(1,:,:))
        case (2)
          Dcoefficients(1,:,:) = -dalpha*SYS_PI*cos(SYS_PI*Dpoints(2,:,:))
        end select
        
      case (3)
        select case(icomp)
        case (1)
          Dcoefficients(1,:,:) = -dalpha*3.0_DP*Dpoints(1,:,:)**2
        case (2)
          Dcoefficients(1,:,:) = -dalpha*3.0_DP*Dpoints(2,:,:)**2
        end select
        
      case (4)
        select case(icomp)
        case (1)
          Dcoefficients(1,:,:) = funcStaticBubble(Dpoints(1,:,:), Dpoints(2,:,:), dalpha, dbeta, dgamma)
        case (2)
          Dcoefficients(1,:,:) = funcStaticBubble(Dpoints(2,:,:), Dpoints(1,:,:), dalpha, dbeta, dgamma)
        end select
      end select
    end select
    
  contains

    elemental real(DP) function funcStaticBubble(dx, dy, dsigma, dh, dr0)
    real(DP), intent(in) :: dx, dy, dsigma, dh, dr0
    
    real(DP), parameter :: dx0 = 0.5_DP
    real(DP), parameter :: dy0 = 0.5_DP
    real(DP) :: dr, dn, dkappa, ddist, dw, ddelta
  
      dr = sqrt((dx-dx0)**2 + (dy-dy0)**2)
      dkappa = -1.0_DP / dr
      dn = (dx-dx0) / dr
      ddist = dr - dr0
      dw = ddist / dh
      if(abs(dw) < 1.0_DP) then
        ddelta = 35.0_DP/32.0_DP * (1.0_DP - 3.0_DP*dw**2 + 3_DP*dw**4 - dw**6) / dh
        funcStaticBubble = ddelta*dkappa*dn*dsigma
      else
        funcStaticBubble = 0.0_DP
      end if

    end function


  end subroutine
  ! ***********************************************************************************************

  subroutine stdrv_funcVelocity2D (icomponent, cderivative, rdiscretisation, &
                nelements, npointsPerElement, Dpoints, rdomainIntSubset, &
                Dvalues, rcollection)
  integer, intent(in)                          :: icomponent
  integer, intent(in)                          :: cderivative
  type(t_spatialDiscretisation), intent(in)    :: rdiscretisation
  integer, intent(in)                          :: nelements
  integer, intent(in)                          :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)       :: Dpoints
  type(t_domainIntSubset), intent(in)          :: rdomainIntSubset
  type(t_collection), intent(inout), optional  :: rcollection
  real(DP), dimension(:,:), intent(out)        :: Dvalues

  integer :: idriver, icomp, ivelo
  real(DP) :: dnu,dc,pi2

    ! fetch coefficients from collection
    dnu = rcollection%DquickAccess(1)
    dc = rcollection%DquickAccess(2)
    idriver = rcollection%IquickAccess(1)
    icomp = rcollection%IquickAccess(5)
    ivelo = rcollection%IquickAccess(6)
    
    pi2 = 2.0_DP * SYS_PI
    
    select case(idriver)
    case (2001)
      select case(icomponent)
      case (1)
        ! X-Velocity
        select case(cderivative)
        case (DER_FUNC2D)
          ! u_1(x,y) = y*(1-y)
          Dvalues(:,:) = Dpoints(2,:,:) * (1.0_DP - Dpoints(2,:,:))
        case (DER_DERIV2D_X)
          ! d_x u_1(x,y) = 0
          Dvalues(:,:) = 0.0_DP
        case (DER_DERIV2D_Y)
          ! d_y u_2(x,y) = 1 - 2*y
          Dvalues(:,:) = 1.0_DP - 2.0_DP * Dpoints(2,:,:)
        end select
      case (2)
        ! Y-Velocity
        ! u_2(x,y) = d_x u_2(x,y) = d_y u_2(x,y) = 0
        Dvalues(:,:) = 0.0_DP
      end select

    case (2002)
      select case(ivelo)
      case (0)
        ! no-flow
        Dvalues(:,:) = 0.0_DP
        
      case (1)
        select case(icomponent)
        case (1)
          ! X-velocity
          select case(cderivative)
          case (DER_FUNC2D)
            ! u_1(x,y) = sin(pi*x) * cos(pi*y)
            Dvalues(:,:) = sin(SYS_PI * Dpoints(1,:,:)) * cos(SYS_PI * Dpoints(2,:,:))
          case (DER_DERIV2D_X)
            ! d_x u_1(x,y) = pi * cos(pi*x) * cos(pi*y)
            Dvalues(:,:) = SYS_PI * cos(SYS_PI * Dpoints(1,:,:)) * cos(SYS_PI * Dpoints(2,:,:))
          case (DER_DERIV2D_Y)
            ! d_y u_1(x,y) = -pi * sin(pi*x) * sin(pi*y)
            Dvalues(:,:) = -SYS_PI * sin(SYS_PI * Dpoints(1,:,:)) * sin(SYS_PI * Dpoints(2,:,:))
          end select
        case (2)
          ! Y-velocity
          select case(cderivative)
          case (DER_FUNC2D)
            ! u_2(x,y) = -cos(pi*x) * sin(pi*y)
            Dvalues(:,:) = -cos(SYS_PI * Dpoints(1,:,:)) * sin(SYS_PI * Dpoints(2,:,:))
          case (DER_DERIV2D_X)
            ! d_x u_2(x,y) = pi * sin(pi*x) * sin(pi*y)
            Dvalues(:,:) = SYS_PI * sin(SYS_PI * Dpoints(1,:,:)) * sin(SYS_PI * Dpoints(2,:,:))
          case (DER_DERIV2D_Y)
            ! d_y u_2(x,y) = -pi * cos(pi*x) * cos(pi*y)
            Dvalues(:,:) = -SYS_PI * cos(SYS_PI * Dpoints(1,:,:)) * cos(SYS_PI * Dpoints(2,:,:))
          end select
        end select
        
      case (2)
        select case(icomponent)
        case (1)
          ! X-velocity
          select case(cderivative)
          case (DER_FUNC2D)
            ! u_1(x,y) = (1 - cos(2*pi*x)) * sin(2*pi*y)
            Dvalues(:,:) = (1.0_DP - cos(pi2 * Dpoints(1,:,:))) * sin(pi2 * Dpoints(2,:,:))
          case (DER_DERIV2D_X)
            ! d_x u_1(x,y) = 2 * pi * sin(2*pi*x) * sin(2*pi*y)
            Dvalues(:,:) = pi2 * sin(pi2 * Dpoints(1,:,:)) * sin(pi2 * Dpoints(2,:,:))
          case (DER_DERIV2D_Y)
            ! d_y u_1(x,y) = 2 * pi *(1- cos(2*pi*x)) * cos(2*pi*y)
            Dvalues(:,:) = pi2 * (1.0_DP - cos(pi2 * Dpoints(1,:,:))) * cos(pi2 * Dpoints(2,:,:))
          end select
        case (2)
          ! Y-velocity
          select case(cderivative)
          case (DER_FUNC2D)
            ! u_2(x,y) = -(1 - cos(2*pi*y)) * sin(2*pi*x)
            Dvalues(:,:) = - (1.0_DP - cos(pi2 * Dpoints(2,:,:))) * sin(pi2 * Dpoints(1,:,:))
          case (DER_DERIV2D_X)
            ! d_x u_2(x,y) = -2 * pi * (1 - cos(pi*y)) * sin(pi*y)
            Dvalues(:,:) = -pi2 * (1.0_DP - cos(pi2 * Dpoints(2,:,:))) * cos(pi2 * Dpoints(1,:,:))
          case (DER_DERIV2D_Y)
            ! d_y u_2(x,y) = -2 * pi * sin(pi*x) * sin(pi*y)
            Dvalues(:,:) = -pi2 * sin(pi2 * Dpoints(1,:,:)) * sin(pi2 * Dpoints(2,:,:))
          end select
        end select
      
      case (3)
        select case(icomponent)
        case (1)
          ! X-Velocity
          select case(cderivative)
          case (DER_FUNC2D)
            Dvalues(:,:) = +2*Dpoints(1,:,:)**2*(1-Dpoints(1,:,:))**2*(Dpoints(2,:,:)*&
                (1-Dpoints(2,:,:))**2-Dpoints(2,:,:)**2*(1-Dpoints(2,:,:)))
          case (DER_DERIV2D_X)
            Dvalues(:,:) = 4*Dpoints(1,:,:)*Dpoints(2,:,:)*(2*Dpoints(2,:,:)-1)*(Dpoints(2,:,:)-1)*&
                (2*Dpoints(1,:,:)-1)*(Dpoints(1,:,:)-1)
          case (DER_DERIV2D_Y)
            Dvalues(:,:) = 2*Dpoints(1,:,:)**2*(Dpoints(1,:,:)-1)**2*(1-6*Dpoints(2,:,:)+6*Dpoints(2,:,:)**2)
          end select
        case (2)
          ! Y-Velocity
          select case(cderivative)
          case (DER_FUNC2D)
            Dvalues(:,:) = -2*Dpoints(2,:,:)**2*(1-Dpoints(2,:,:))**2*(Dpoints(1,:,:)*&
                (1-Dpoints(1,:,:))**2-Dpoints(1,:,:)**2*(1-Dpoints(1,:,:)))
          case (DER_DERIV2D_X)
            Dvalues(:,:) = -2*Dpoints(2,:,:)**2*(Dpoints(2,:,:)-1)**2*(1-6*Dpoints(1,:,:)+6*Dpoints(1,:,:)**2)
          case (DER_DERIV2D_Y)
            Dvalues(:,:) = -4*Dpoints(1,:,:)*Dpoints(2,:,:)*(2*Dpoints(2,:,:)-1)*(Dpoints(2,:,:)-1)*&
                (2*Dpoints(1,:,:)-1)*(Dpoints(1,:,:)-1)
          end select
        end select
      end select
    end select

  end subroutine

  ! ***********************************************************************************************

  subroutine stdrv_funcPressure2D (icomponent, cderivative, rdiscretisation, &
                nelements, npointsPerElement, Dpoints, rdomainIntSubset, &
                Dvalues, rcollection)
  integer, intent(in)                          :: icomponent
  integer, intent(in)                          :: cderivative
  type(t_spatialDiscretisation), intent(in)    :: rdiscretisation
  integer, intent(in)                          :: nelements
  integer, intent(in)                          :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)       :: Dpoints
  type(t_domainIntSubset), intent(in)          :: rdomainIntSubset
  type(t_collection), intent(inout), optional  :: rcollection
  real(DP), dimension(:,:), intent(out)        :: Dvalues

  integer :: idriver, icomp, ipres
  real(DP) :: dnu,dalpha,dpi2

    ! fetch coefficients from collection
    dnu = rcollection%DquickAccess(1)
    dalpha = rcollection%DquickAccess(2)
    idriver = rcollection%IquickAccess(1)
    icomp = rcollection%IquickAccess(5)
    ipres = rcollection%IquickAccess(7)
    
    select case(idriver)
    case (2001)
      ! p(x,y) = 2*nu*alpha*(1-x)
      Dvalues(:,:) = 2.0_DP*dnu*dalpha*(1.0_DP - Dpoints(1,:,:))
      
    case (2002)
      select case(ipres)
      case (0)
        ! p(x,y) = 0
        Dvalues(:,:) = 0.0_DP
        
      case (1)
        ! p(x,y) = alpha * cos(pi*x) * cos(pi*y)
        Dvalues(:,:) = dalpha * cos(SYS_PI * Dpoints(1,:,:)) * cos(SYS_PI * Dpoints(2,:,:))
        
      case (2)
        ! p(x,y) = alpha*(4/pi - sin(pi*x) - sin(pi*y))
        Dvalues(:,:) = dalpha*(4.0_DP / SYS_PI - sin(SYS_PI * Dpoints(1,:,:)) - sin(SYS_PI * Dpoints(2,:,:)))
      
      case (3)
        ! p(x,y) = alpha*(1/2 - x^3 - y^3)
        Dvalues(:,:) = dalpha*(0.5_DP - Dpoints(1,:,:)**3 - Dpoints(2,:,:)**3)

      end select
      
    end select
    
  end subroutine
  
end module
