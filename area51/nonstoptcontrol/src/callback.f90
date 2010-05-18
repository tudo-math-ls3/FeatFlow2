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
  
    rphysics%cequation = 0
    rphysics%dviscosity = 1.0_DP
    rphysics%doptControlAlpha = 1.0_DP
    rphysics%doptControlGamma = 0.0_DP
    rphysics%dcouplePrimalToDual = 1.0_DP
    rphysics%dcoupleDualToPrimal = 1.0_DP
  
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
    real(DP) :: dtime,dalpha,dtstep
    
    icomponent = rcollection%IquickAccess(1)
    dtime = rcollection%DquickAccess(1)
    dalpha = rcollection%DquickAccess(2)
    dtstep = rcollection%DquickAccess(3)
    
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
        !                     + 2*dtime**2*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))) 
        
        !                     + ((1.0_DP-dtime)**2/dalpha)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
                            
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
      Dcoefficients(1,:,:) = 2.0_DP*Dpoints(1,:,:) - dtime**2*Dpoints(1,:,:) 
      
      if (dtime .eq. 1.0_DP) then
        Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
      end if
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
    
    call boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    
    Dvalues(:) = 0.0_DP

    select case (cequation)
    case (0)
      ! Heat equation
      select case (icomponent)
      case (1)
        !Dvalues(1) = dtime*(dx*dy)
        !Dvalues(1) = 0.0_DP
        Dvalues(1) = dtime**2 * dx
      case (2)
        !Dvalues(1) = (1.0_DP-dtime)*(dx*dy)
        !Dvalues(1) = 0.0_DP
        !Dvalues(1) = - (2.0_DP*dtime**2 * dy*(1.0_DP-dy)  +  2.0_DP*dtime**2.0_DP * dx*(1.0_DP-dx) )
        Dvalues(1) = -2.0_DP*dtime*dx
      case default
        ! Should not happen
        call sys_halt()
      end select
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
