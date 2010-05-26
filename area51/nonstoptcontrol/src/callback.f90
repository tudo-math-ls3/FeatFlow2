!##############################################################################
!# ****************************************************************************
!# <name> callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# Callback routines. Boundary conditions, etc...
!# </purpose>
!##############################################################################

module callback

  use fsystem
  use storage
  use genoutput
  use basicgeometry
  use collection
  use spatialdiscretisation
  use boundary
  use discretebc
  
  use scalarpde
  use linearformevaluation
  use collection
  use domainintegration
  
  use physics
  
  implicit none

contains

  ! ***************************************************************************

  subroutine cb_getPhysicsHeatEqn(rphysics)
  
  ! Initialises the physics structure for the equation
  
  ! The physics structure to initialise
  type(t_physics), intent(out) :: rphysics
  
    ! 0 = 2D heat equation
    ! 1 = 2D Stokes
    ! 2 = 1D heat equation
    rphysics%cequation = 2
    rphysics%dviscosity = 1.0_DP
    rphysics%doptControlAlpha = 1.0_DP
    rphysics%doptControlGamma = 0.0_DP
    rphysics%dcouplePrimalToDual = 1.0_DP
    rphysics%dcoupleDualToPrimal = 1.0_DP
    rphysics%dcoupleTermCond = 1.0_DP
  
  end subroutine

  ! ***************************************************************************

  subroutine rhs_heatEquation (rdiscretisation, rform, &
                nelements, npointsPerElement, Dpoints, &
                IdofsTest, rdomainIntSubset, &
                Dcoefficients, rcollection)
  
  ! Coefficient function for the RHS of the heat equation
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  type(t_linearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset    
  type(t_collection), intent(inout), optional :: rcollection

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  
    integer :: icomponent
    real(DP) :: dtime,dalpha,dtstep,dgamma
    
    icomponent = rcollection%IquickAccess(1)
    dtime = rcollection%DquickAccess(1)
    dtstep = rcollection%DquickAccess(2)
    dalpha = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    
    Dcoefficients(1,:,:) = 0.0_DP
    
    ! Primal or dual?
    select case (icomponent)
    case (1)
      ! Primal: f
      if (dtime .gt. 0.0_DP) then
        !Dcoefficients(1,:,:) = Dpoints(1,:,:)*Dpoints(2,:,:) &
        !                     + ((1.0_DP-dtime)/dalpha)*Dpoints(1,:,:)*Dpoints(2,:,:)
        !Dcoefficients(1,:,:) = Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
        !                     + 2*dtime*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))) 
        !                     + ((1.0_DP-dtime)/dalpha)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
        
        !Dcoefficients(1,:,:) = 2*dtime*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
        !                     + 2*dtime**2*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)))        !                     + ((1.0_DP-dtime)**2/dalpha)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
        
        ! 1.)
        ! Dcoefficients(1,:,:) = 0.0_DP
        
        ! 2.)
        !Dcoefficients(1,:,:) = &
        !  3.0_DP*dtime*Dpoints(1,:,:)-7.0_DP*dtime**2*Dpoints(1,:,:)+4*dtime**3*Dpoints(1,:,:)

        ! 3.)
        Dcoefficients(1,:,:) = &
          3.0_DP*dtime*Dpoints(1,:,:)-8.0_DP*dtime**2*Dpoints(1,:,:)+5*dtime**3*Dpoints(1,:,:)
                            
      end if
      
    case (2)
      ! Dual: -z
      !Dcoefficients(1,:,:) = Dpoints(1,:,:)*Dpoints(2,:,:) &
      !                     -dtime*Dpoints(1,:,:)*Dpoints(2,:,:)
      !Dcoefficients(1,:,:) = Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
      !                     + 2*(1.0_DP-dtime)*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))) 
      !                      - dtime*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
      !Dcoefficients(1,:,:) = 2*(dtime-1.0_DP)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
      !                     + 2*(1.0_DP-dtime)**2*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))) 
      !                       - dtime**2*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
      !Dcoefficients(1,:,:) = (2.0_DP*Dpoints(1,:,:)*(1-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
      !                     - 8.0_DP*dtime**2 &
      !                     + dtime**2*Dpoints(1,:,:)*(1-Dpoints(1,:,:))*Dpoints(2,:,:)*(1-Dpoints(2,:,:)) )
      
      !Dcoefficients(1,:,:) = -(-2.0_DP*Dpoints(1,:,:) + (dtime**2-2*dtime)*Dpoints(1,:,:))
      
      ! 1.)
      !Dcoefficients(1,:,:) = -(-2.0_DP*Dpoints(1,:,:) + dtime**2*Dpoints(1,:,:))
      
      ! 2.)
      Dcoefficients(1,:,:) = &
          -( (1.0_DP-dtime)**2*Dpoints(1,:,:)-2.0_DP*dtime*(1.0_DP-dtime)*Dpoints(1,:,:)+&
             dtime**2*(1.0_DP-dtime)**2*Dpoints(1,:,:) )
      
      if (dtime .eq. 1.0_DP) then
        !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
        Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP + dgamma/dtstep)
      end if
    case default
      ! Should not happen
      call sys_halt()
    end select      

  end subroutine
  
  ! ***************************************************************************

  subroutine rhs_stokesEquation (rdiscretisation, rform, &
                nelements, npointsPerElement, Dpoints, &
                IdofsTest, rdomainIntSubset, &
                Dcoefficients, rcollection)
  
  ! Coefficient function for the RHS of the heat equation
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  type(t_linearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset    
  type(t_collection), intent(inout), optional :: rcollection

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  
    integer :: icomponent
    real(DP) :: dtime,dalpha,dtstep
    
    icomponent = rcollection%IquickAccess(1)
    dtime = rcollection%DquickAccess(1)
    dalpha = rcollection%DquickAccess(2)
    dtstep = rcollection%DquickAccess(3)
    
    Dcoefficients(1,:,:) = 0.0_DP
    
    ! Primal or dual?
    select case (icomponent)
    case (1)
      ! Primal: f1

    case (2)
      ! Primal: f2
      
    case (4)
      ! Dual: -z1
      
      ! 1.)
      Dcoefficients(1,:,:) =  2.0_DP*Dpoints(1,:,:) - dtime**2*Dpoints(1,:,:) 

    case (5)
      ! Dual: -z2
      
      ! 1.)
      Dcoefficients(1,:,:) =  - 2.0_DP*Dpoints(2,:,:) + dtime**2*Dpoints(2,:,:) 

    case default
      ! Should not happen
      call sys_halt()
    end select      


  end subroutine
  
  ! ***************************************************************************
  
  subroutine cb_getBoundaryValuesOptC (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
      cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
                                  
  ! Calculates values of the primal solution on the boundary
  
  integer, dimension(:), intent(in)                           :: Icomponents
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  integer, intent(in)                                         :: ielement
  integer, intent(in)                                         :: cinfoNeeded
  integer, intent(in)                                         :: iwhere
  real(DP), intent(in)                                        :: dwhere
  type(t_collection), intent(inout), optional                 :: rcollection

  real(DP), dimension(:), intent(out)                         :: Dvalues

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    real(DP) :: dx,dy
    real(DP) :: dtime
    integer :: icomponent,cequation
    
    icomponent = rcollection%IquickAccess(1)
    cequation = rcollection%IquickAccess(2)
    dtime = rcollection%DquickAccess(1)
    
    if (rdiscretisation%p_rtriangulation%ndim .eq. NDIM2D) then
      call boundary_getCoords(rdiscretisation%p_rboundary, &
          rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    else
      dx = dwhere
      dy = 0.0_DP
    end if
    
    Dvalues(:) = 0.0_DP

    select case (cequation)
    case (0,2)
      ! Heat equation
      select case (icomponent)
      case (1)
        !Dvalues(1) = dtime*(dx*dy)
        !Dvalues(1) = 0.0_DP
        
        ! 1.)
        !Dvalues(1) = dtime**2 * dx 
        
        ! 2.) -> BC in spacetimebc.f90 beachten!
        !Dvalues(1) = dtime**2 * (1.0_DP-dtime)**2 * dx 

        ! 3.) -> BC in spacetimebc.f90 beachten!
        Dvalues(1) = dtime**2 * (1.0_DP-dtime)**2 * dx 
      case (2)
        !Dvalues(1) = (1.0_DP-dtime)*(dx*dy)
        !Dvalues(1) = 0.0_DP
        !Dvalues(1) = - (2.0_DP*dtime**2 * dy*(1.0_DP-dy)  +  2.0_DP*dtime**2.0_DP * dx*(1.0_DP-dx) )
        
        ! 1.)
        !Dvalues(1) = - 2.0_DP*dtime * dx
        
        ! 2.) -> BC in spacetimebc.f90 beachten!
        !Dvalues(1) = dtime * (1.0_DP-dtime) * dx 

        ! 3.) -> BC in spacetimebc.f90 beachten!
        Dvalues(1) = dtime * (1.0_DP-dtime)**2 * dx 
      case default
        ! Should not happen
        call sys_halt()
      end select

    case (1)
      ! Stokes equation
      select case (icomponent)
      
      ! Primal BC
      case (1)
        Dvalues(1) = dtime**2 * dx
      case (2)
        Dvalues(1) = - dtime**2 * dy
      
      ! Dual BC
      case (4)
      
        ! 1.)
        Dvalues(1) = -2.0_DP*dtime*dx
        
        
      case (5)
      
        ! 1.)
        Dvalues(1) = 2.0_DP*dtime*dy
      
      case default
        ! Should not happen
        call sys_halt()
      end select
      
    case default
    
      call output_line ("Equation not supported.")
      call sys_halt()
      
    end select

  end subroutine


  subroutine ffunctionTermCond (cderivative, rdiscretisation, &
                                  nelements, npointsPerElement, Dpoints, &
                                  IdofsTest, rdomainIntSubset, &
                                  Dvalues, rcollection)
  
  ! Terminal condition
  
  integer, intent(in) :: cderivative
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  type(t_collection), intent(inout), optional :: rcollection
  
  real(DP), dimension(:,:), intent(out) :: Dvalues

    Dvalues(:,:) = -2.0_DP*Dpoints(1,:,:)

  end subroutine

end module
