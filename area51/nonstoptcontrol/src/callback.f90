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
  use paramlist
  
  use physics
  
  implicit none

contains

  ! ***************************************************************************

  subroutine cb_getPhysicsHeatEqn(rparlist,rphysics)
  
  ! Initialises the physics structure for the equation
  
  ! Parameter list containing the DAT file parameters
  type(t_parlist), intent(in) :: rparlist

  ! The physics structure to initialise
  type(t_physics), intent(out) :: rphysics
  
    ! 0 = 2D heat equation
    ! 1 = 2D Stokes
    ! 2 = 1D heat equation
!    rphysics%cequation = 2
!    rphysics%creferenceProblem = 4
!    rphysics%dviscosity = 1.0_DP
!    rphysics%doptControlAlpha = 0.01_DP
!    rphysics%doptControlGamma = 0.0_DP
!    rphysics%dcouplePrimalToDual = 1.0_DP
!    rphysics%dcoupleDualToPrimal = 1.0_DP
!    rphysics%dcoupleTermCond = 1.0_DP
  
    call parlst_getvalue_int (rparlist, "PHYSICS", &
        "cequation", rphysics%cequation)
    call parlst_getvalue_int (rparlist, "PHYSICS", &
        "creferenceProblem", rphysics%creferenceProblem)
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "dviscosity", rphysics%dviscosity)
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "doptControlAlpha", rphysics%doptControlAlpha)
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "doptControlGamma", rphysics%doptControlGamma)
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "dcouplePrimalToDual", rphysics%dcouplePrimalToDual)
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "dcoupleDualToPrimal", rphysics%dcoupleDualToPrimal)
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "dcoupleTermCond", rphysics%dcoupleTermCond)
  
  end subroutine

  ! ***************************************************************************
  ! Reference functions, Heat equation
  ! ***************************************************************************
  ! ***************************************************************************
  ! Heat equation, function set 1.
  !   y = t^2   x1 
  ! ***************************************************************************

  elemental real(DP) function fct_heatY1 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatY1 = dtime**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda1 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatLambda1 = dalpha * (- 2.0_DP*dtime * dx) + 2*dalpha*dx
  end function
    
  elemental real(DP) function fct_heatF1 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatF1 = 2.0_DP*dtime*dx - 2.0_DP*dx*(dtime-1.0_DP)
  end function
  
  elemental real(DP) function fct_heatZ1 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatZ1 = -2.0_DP*dalpha*dx + dtime**2*dx
  end function

  ! ***************************************************************************
  ! Heat equation, function set 2
  !   y = t^2 (1-t)^2   x1 
  !   l = a t (1-t) x1
  ! No factor 4 due to inconsistent terminal condition.
  ! ***************************************************************************
  
  elemental real(DP) function fct_heatY2 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatY2 = dtime**2 * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda2 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatLambda2 = dalpha * dtime * (1.0_DP-dtime) * dx 
  end function

  elemental real(DP) function fct_heatF2 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatF2 = 3.0_DP*dtime*dx-7.0_DP*dtime**2*dx+4*dtime**3*dx
  end function

  elemental real(DP) function fct_heatZ2 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatZ2 = (1.0_DP-dtime)*dx-dtime*dx+dtime**2*(1.0_DP-dtime)**2*dx 
  end function

  ! ***************************************************************************
  ! Heat equation, function set 3
  !   y = t^2 (1-t)^2   x1 
  !   l = a t (1-t)^2   x1
  ! Gives factor 4 with CN due to proper terminal condition.
  ! ***************************************************************************

  elemental real(DP) function fct_heatY3 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatY3 = dtime**2 * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda3 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatLambda3 = dalpha * dtime * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatF3 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatF3 = &
      3.0_DP*dtime*dx-8.0_DP*dtime**2*dx+5*dtime**3*dx
  end function

  elemental real(DP) function fct_heatZ3 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatZ3 = &
         (1.0_DP-dtime)**2*dx-2.0_DP*dtime*(1.0_DP-dtime)*dx+&
          dtime**2*(1.0_DP-dtime)**2*dx 
  end function
          
  ! ***************************************************************************
  ! Heat equation, function set 4
  !   y = t^2 (1-t)^2   x1 
  ! ***************************************************************************
    
  elemental real(DP) function fct_heatY4 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatY4 = dtime**2 * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda4 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatLambda4 = dalpha * (-2.0_DP*dtime*(1.0_DP-dtime)**2*dx + 2.0_DP*dtime**2*(1.0_DP-dtime)*dx)
  end function
    
  elemental real(DP) function fct_heatF4 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatF4 = 0.0_DP
  end function

  elemental real(DP) function fct_heatZ4 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatZ4 = &
        - dalpha*2.0_DP*(1-dtime)**2*dx + dalpha*8.0_DP*dtime*(1-dtime)*dx &
        - dalpha*2.0_DP*dtime**2*dx + dtime**2*(1.0_DP-dtime)**2*dx
  end function

  ! ***************************************************************************
  ! Heat equation, function set 5.
  !   y = t   x1 
  !   l = (1-t) x1
  ! Does not produce error=0, only at the beginning.
  ! ***************************************************************************

  elemental real(DP) function fct_heatY5 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatY5 = dtime * dx 
  end function

  elemental real(DP) function fct_heatLambda5 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatLambda5 = (1.0_DP-dtime)*dx
  end function
    
  elemental real(DP) function fct_heatF5 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatF5 = -(dtime-dalpha-1.0_DP)*dx/dalpha
  end function
  
  elemental real(DP) function fct_heatZ5 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatZ5 = (dtime-1.0_DP)*dx
  end function

  ! ***************************************************************************
  ! Heat equation, function set 6.
  !   y = t^2   x1 
  !   l = (1-t)^2 x1
  ! Factor 4 in CN due to appropriate terminal condition.
  ! ***************************************************************************

  elemental real(DP) function fct_heatY6 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatY6 = dtime**2 * dx
  end function

  elemental real(DP) function fct_heatLambda6 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatLambda6 = (1.0_DP-dtime)**2*dx
  end function
    
  elemental real(DP) function fct_heatF6 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatF6 = dx*(2.0_DP*dtime*dalpha+1.0_DP-2.0_DP*dtime+dtime**2)/dalpha
  end function
  
  elemental real(DP) function fct_heatZ6 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_heatZ6 = dtime**2*dx-2.0_DP*(1.0_DP-dtime)*dx
  end function

  ! ***************************************************************************
  ! Reference functions, Stokes equation
  ! ***************************************************************************
  ! ***************************************************************************
  ! Stokes, function set 1.
  !   y = t^2  ( x1 , -x2 , 0)
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY1_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY1_x = dtime**2 * dx
  end function

  elemental real(DP) function fct_stokesY1_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY1_y = - dtime**2 * dy
  end function

  elemental real(DP) function fct_stokesP1 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesP1 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda1_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda1_x = dalpha * (-2.0_DP*dtime*dx)
  end function

  elemental real(DP) function fct_stokesLambda1_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda1_y = dalpha * (2.0_DP*dtime*dy)
  end function

  elemental real(DP) function fct_stokesXi1 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesXi1 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ1_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ1_x = -2.0_DP*dx + dtime**2*dx
  end function

  elemental real(DP) function fct_stokesZ1_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ1_y = 2.0_DP*dy - dtime**2*dy
  end function

  ! ***************************************************************************
  ! Stokes, function set 2.
  !   y = t^2 (1-t)^3  ( x2(1-x2) , 0 , 0 )
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY2_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY2_x = dy * (-1 + dy) * dtime ** 2 * (-1 + dtime) ** 3
  end function

  elemental real(DP) function fct_stokesY2_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY2_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP2 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesP2 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda2_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda2_x = &
        - dalpha * dtime * (-1 + dtime) ** 2 * (2 * dtime - 2 * dtime ** 2 &
        + 2 * dy - 5 * dy * dtime - 2 * dy ** 2 + 5 * dy ** 2 * dtime)
  end function

  elemental real(DP) function fct_stokesLambda2_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda2_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi2 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesXi2 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ2_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ2_x = &
          dy * dtime ** 2 - 3 * dy * dtime ** 3 + 3 * dy * dtime ** 4 &
        - dy * dtime ** 5 - dy ** 2 * dtime ** 2 + 3 * dy ** 2 * dtime ** 3 &
        - 3 * dy ** 2 * dtime ** 4 + dy ** 2 * dtime ** 5 - 2 * dalpha * dy &
        - 0.18D2 * dalpha * dy ** 2 * dtime + 2 * dalpha * dy ** 2 &
        + 0.18D2 * dalpha * dy * dtime - 0.36D2 * dalpha * dy * dtime ** 2 &
        + 0.20D2 * dalpha * dy * dtime ** 3 + 0.36D2 * dalpha * dy ** 2 * dtime ** 2 &
        - 0.20D2 * dalpha * dy ** 2 * dtime ** 3
  end function

  elemental real(DP) function fct_stokesZ2_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ2_y = 0.0_DP
  end function

  ! ***************************************************************************
  ! Stokes, function set 3.
  !   y = t^2 (x1^2 x2 , -x1 x2^2 , 0)
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY3_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY3_x = dtime**2 * dx**2 * dy
  end function

  elemental real(DP) function fct_stokesY3_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY3_y = - dtime**2 * dx**2 * dy
  end function

  elemental real(DP) function fct_stokesP3 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesP3 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda3_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda3_x = -dalpha*(-2*dy*dtime**2*(1-dtime)**2 &
        +2*dx**2*dy*dtime*(1-dtime)**2 &
        -2*dx**2*dy*dtime**2*(1-dtime) )
  end function

  elemental real(DP) function fct_stokesLambda3_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda3_y = -dalpha*(-2*dx*dtime**2*(1-dtime)**2 &
        +2*dx*dy**2*dtime*(1-dtime)**2 &
        -2*dx*dy**2*dtime**2*(1-dtime) )
  end function

  elemental real(DP) function fct_stokesXi3 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesXi3 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ3_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ3_x = (dx**2*dy*dtime**2*(-1+dtime)**2) &
        - dalpha*(-4*dy*dtime+12*dy*dtime**2 &
                  -8*dy*dtime**3+2*dx**2*dy &
                  -12*dx**2*dy*dtime+12*dx**2*dy*dtime**2)
  end function

  elemental real(DP) function fct_stokesZ3_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ3_y = (dx*dy**2*dtime**2*(-1+dtime)**2) &
          - dalpha*(-4*dx*dtime+12*dx*dtime**2 &
                    -8*dx*dtime**3+2*dx*dy**2 &
                    -12*dx*dy**2*dtime+12*dx*dy**2*dtime**2)
  end function

  ! ***************************************************************************
  ! Stokes, function set 4.
  !   y = t^2 (1-t)^2  (x1^3 x2^2 , -x1^2 x2^3 , 0)
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY4_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY4_x = dx ** 3 * dy ** 2 * dtime ** 2 * (-1.0_DP + dtime) ** 2
  end function

  elemental real(DP) function fct_stokesY4_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY4_y = -dx ** 2 * dy ** 3 * dtime ** 2 * (-1.0_DP + dtime) ** 2
  end function

  elemental real(DP) function fct_stokesP4 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesP4 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda4_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda4_x = -2 * dalpha * dx * dtime * (-1 + dtime) * &
              (3 * dy ** 2 * dtime - 3 * dy ** 2 * dtime ** 2 + dx ** 2 * dtime &
              - dx ** 2 * dtime ** 2 - dx ** 2 * dy ** 2 + 2 * dx ** 2 * dy ** 2 * dtime)
  end function

  elemental real(DP) function fct_stokesLambda4_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda4_y = 2 * dalpha * dy * dtime * (-1 + dtime) * &
              (dy ** 2 * dtime - dy ** 2 * dtime ** 2 + 3 * dx ** 2 * dtime &
               - 3 * dx ** 2 * dtime ** 2 - dx ** 2 * dy ** 2 &
               + 2 * dx ** 2 * dy ** 2 * dtime)
  end function

  elemental real(DP) function fct_stokesXi4 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesXi4 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ4_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ4_x = 6 * dx * dy ** 2 * dtime ** 2 &
        - 0.12D2 * dx * dy ** 2 * dtime ** 3 &
        + 6 * dx * dy ** 2 * dtime ** 4 + 2 * dx ** 3 * dtime ** 2 &
        - 4 * dx ** 3 * dtime ** 3 + 2 * dx ** 3 * dtime ** 4 &
        + dx ** 3 * dy ** 2 * dtime ** 2 - 2 * dx ** 3 * dy ** 2 * dtime ** 3 &
        + dx ** 3 * dy ** 2 * dtime ** 4 + 0.12D2 * dalpha * dx * dy ** 2 * dtime &
        - 0.36D2 * dalpha * dx * dy ** 2 * dtime ** 2 &
        + 0.24D2 * dalpha * dx * dy ** 2 * dtime ** 3 &
        + 4 * dalpha * dx ** 3 * dtime &
        - 0.12D2 * dalpha * dx ** 3 * dtime ** 2 &
        + 8 * dalpha * dx ** 3 * dtime ** 3 &
        - 2 * dalpha * dx ** 3 * dy ** 2 &
        + 0.12D2 * dalpha * dx ** 3 * dy ** 2 * dtime &
        - 0.12D2 * dalpha * dx ** 3 * dy ** 2 * dtime ** 2
  end function

  elemental real(DP) function fct_stokesZ4_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ4_y = -2 * dy ** 3 * dtime ** 2 + 4 * dy ** 3 * dtime ** 3 &
        - 2 * dy ** 3 * dtime ** 4 - 6 * dx ** 2 * dy * dtime ** 2 &
        + 0.12D2 * dx ** 2 * dy * dtime ** 3 - 6 * dx ** 2 * dy * dtime ** 4 &
        - dx ** 2 * dy ** 3 * dtime ** 2 + 2 * dx ** 2 * dy ** 3 * dtime ** 3 &
        - dx ** 2 * dy ** 3 * dtime ** 4 - 4 * dalpha * dy ** 3 * dtime &
        + 0.12D2 * dalpha * dy ** 3 * dtime ** 2 &
        - 8 * dalpha * dy ** 3 * dtime ** 3 &
        - 0.12D2 * dalpha * dx ** 2 * dy * dtime &
        + 0.36D2 * dalpha * dx ** 2 * dy * dtime ** 2 &
        - 0.24D2 * dalpha * dx ** 2 * dy * dtime ** 3 &
        + 2 * dalpha * dx ** 2 * dy ** 3 &
        - 0.12D2 * dalpha * dx ** 2 * dy ** 3 * dtime &
        + 0.12D2 * dalpha * dx ** 2 * dy ** 3 * dtime ** 2
  end function
  
  ! ***************************************************************************
  ! Stokes, function set 5.
  !   y = t  ( x2(1-x2) , 0 , 0 )
  ! No coupling between primal and dual!
  ! EXAMPLE DOES NOT WORK, NOT DIVERGENCE FREE WHEN BEING TIME DEPENDENT!
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY5_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY5_x = dtime * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesY5_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY5_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP5 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesP5 = dtime*(1.0_DP-dx) !0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda5_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda5_x = 0.0_DP !-dalpha*( &
!                2.0_DP*dtime**2*(1.0_DP-dtime)**2 &
!              + 2.0_DP*dy*(1.0_DP-dy)*dtime*(1.0_DP-dtime)**2 &
!              - 2.0_DP*dy*(1.0_DP-dy)*dtime**2*(1.0_DP-dtime))
  end function

  elemental real(DP) function fct_stokesLambda5_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda5_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi5 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesXi5 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ5_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ5_x = 0.0_DP !- dy*(-1.0_DP+dy)*dtime**2*(-1.0_DP+dtime)**2 &
!            - dalpha*(4.0_DP*dtime-12.0_DP*dtime**2+8.0_DP*dtime**3+2.0_DP*dy &
!                      -12.0_DP*dy*dtime+12.0_DP*dy*dtime**2 &
!                      -2.0_DP*dy**2+12.0_DP*dy**2*dtime &
!                      -12.0_DP*dy**2*dtime**2)
  end function

  elemental real(DP) function fct_stokesZ5_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ5_y = 0.0_DP
  end function

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
  
    integer :: icomponent,creferenceProblem
    real(DP) :: dtime,dalpha,dtstep,dgamma,dtheta
    
    icomponent = rcollection%IquickAccess(1)
    creferenceProblem = rcollection%IquickAccess(2)
    dtime = rcollection%DquickAccess(1)
    dtstep = rcollection%DquickAccess(2)
    dalpha = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    dtheta = rcollection%DquickAccess(5)
    
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
        
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dcoefficients(1,:,:) = fct_heatF1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) 
        
        case (2)
          ! 2.)
          Dcoefficients(1,:,:) = fct_heatF2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) 

        case (3)
          ! 3.)
          Dcoefficients(1,:,:) = fct_heatF3(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) 

        case (4)
          ! 4.)
          Dcoefficients(1,:,:) = fct_heatF4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) 

        case (5)
          ! 5.)
          Dcoefficients(1,:,:) = fct_heatF5(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) 

        case (6)
          ! 6.)
          Dcoefficients(1,:,:) = fct_heatF6(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) 

        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
                            
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
      
      select case (creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        Dcoefficients(1,:,:) = &
            - (fct_heatZ1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      
      case (2)  
        ! 2.)
        Dcoefficients(1,:,:) = &
            - (fct_heatZ2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )

      case (3)  
        ! 3.)
        Dcoefficients(1,:,:) = &
            - (fct_heatZ3(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
              
      case (4)
        ! 4.)
        Dcoefficients(1,:,:) = &
            - (fct_heatZ4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )

      case (5)
        ! 5.)
        Dcoefficients(1,:,:) = &
            - (fct_heatZ5(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )

      case (6)
        ! 6.)
        Dcoefficients(1,:,:) = &
            - (fct_heatZ6(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )

      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select
      
      ! Dual solution lives in the interval midpoints.
      ! Dual RHS lives in the endpoints of the interval!
      ! dtime = 1.0 is ok!
      if (dtime .eq. 1.0_DP) then
        !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
        Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (dtheta + dgamma/dtstep)
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
  
    integer :: icomponent,creferenceProblem
    real(DP) :: dtime,dalpha,dtstep,dtheta,dgamma
    
    icomponent = rcollection%IquickAccess(1)
    creferenceProblem = rcollection%IquickAccess(2)
    dtime = rcollection%DquickAccess(1)
    dtstep = rcollection%DquickAccess(2)
    dalpha = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    dtheta = rcollection%DquickAccess(5)
    
    Dcoefficients(1,:,:) = 0.0_DP
    
    ! Primal or dual?
    select case (icomponent)
    case (1)
      ! Primal: f1

    case (2)
      ! Primal: f2
      
    case (4)
      ! Dual: -z1
      
      select case (creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        Dcoefficients(1,:,:) =  &
            - (fct_stokesZ1_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      case (2)
        ! 2.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ2_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      case (3)
        ! 3.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ3_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      
      case (4)
        ! 4.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )

      case (5)
        ! 5.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ5_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      
      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select

      if (dtime .eq. 1.0_DP) then
        !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
        Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (dtheta + dgamma/dtstep)
      end if

    case (5)
      ! Dual: -z2
      
      select case (creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        Dcoefficients(1,:,:) =  &
            - (fct_stokesZ1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      case (2)
        ! 2.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      case (3)
        ! 3.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ3_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      case (4)
        ! 4.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      case (5)
        ! 5.)
        Dcoefficients(1,:,:) = &
            - (fct_stokesZ5_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha) )
      
      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select

      if (dtime .eq. 1.0_DP) then
        !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
        Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (dtheta + dgamma/dtstep)
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
    real(DP) :: dtime,dalpha
    integer :: icomponent,cequation,creferenceProblem
    
    icomponent = rcollection%IquickAccess(1)
    cequation = rcollection%IquickAccess(2)
    creferenceProblem = rcollection%IquickAccess(3)
    dtime = rcollection%DquickAccess(1)
    dalpha = rcollection%DquickAccess(2)
    
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
      ! -------------------------------------------------------------
      ! Heat equation
      ! -------------------------------------------------------------
      select case (icomponent)
      case (1)
        !Dvalues(1) = dtime*(dx*dy)
        !Dvalues(1) = 0.0_DP
        
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_heatY1 (dx,dy,dtime,dalpha)
        
        case (2)
          ! 2.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY2 (dx,dy,dtime,dalpha)

        case (3)
          ! 3.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY3 (dx,dy,dtime,dalpha)

        case (4)
          ! 4.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY4 (dx,dy,dtime,dalpha)

        case (5)
          ! 5.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY5 (dx,dy,dtime,dalpha)

        case (6)
          ! 6.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY6 (dx,dy,dtime,dalpha)

        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
        
      case (2)
        !Dvalues(1) = (1.0_DP-dtime)*(dx*dy)
        !Dvalues(1) = 0.0_DP
        !Dvalues(1) = - (2.0_DP*dtime**2 * dy*(1.0_DP-dy)  +  2.0_DP*dtime**2.0_DP * dx*(1.0_DP-dx) )
        
        ! For 1D: See BC in spacetimebc.f90!!!
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_heatLambda1 (dx,dy,dtime,dalpha)
          
        case (2)
          ! 2.) 
          Dvalues(1) = fct_heatLambda2 (dx,dy,dtime,dalpha)

        case (3)
          ! 3.) 
          Dvalues(1) = fct_heatLambda3 (dx,dy,dtime,dalpha)
          
        case (4)
          ! 4.) 
          Dvalues(1) = fct_heatLambda4 (dx,dy,dtime,dalpha)
          
        case (5)
          ! 5.) 
          Dvalues(1) = fct_heatLambda5 (dx,dy,dtime,dalpha)
          
        case (6)
          ! 6.) 
          Dvalues(1) = fct_heatLambda6 (dx,dy,dtime,dalpha)
          
        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
      case default
        ! Should not happen
        call sys_halt()
      end select

    case (1)
      ! -------------------------------------------------------------
      ! Stokes equation
      ! -------------------------------------------------------------
      select case (icomponent)
      
      ! Primal BC
      case (1)
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_stokesY1_x (dx,dy,dtime,dalpha)
        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesY2_x (dx,dy,dtime,dalpha)
        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesY3_x (dx,dy,dtime,dalpha)

        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesY4_x (dx,dy,dtime,dalpha)

        case (5)
          ! 5.)
          Dvalues(1) = fct_stokesY5_x (dx,dy,dtime,dalpha)

        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
        
      case (2)
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_stokesY1_y (dx,dy,dtime,dalpha)

        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesY2_y (dx,dy,dtime,dalpha)

        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesY3_y (dx,dy,dtime,dalpha)

        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesY4_y (dx,dy,dtime,dalpha)

        case (5)
          ! 4.)
          Dvalues(1) = fct_stokesY5_y (dx,dy,dtime,dalpha)

        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
      
      ! Dual BC
      case (4)
      
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_stokesLambda1_x (dx,dy,dtime,dalpha)
        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesLambda2_x (dx,dy,dtime,dalpha)
        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesLambda3_x (dx,dy,dtime,dalpha)

        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesLambda4_x (dx,dy,dtime,dalpha)

        case (5)
          ! 5.)
          Dvalues(1) = fct_stokesLambda5_x (dx,dy,dtime,dalpha)
        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
        
      case (5)
      
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_stokesLambda1_y (dx,dy,dtime,dalpha)
        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesLambda2_y (dx,dy,dtime,dalpha)
        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesLambda3_y (dx,dy,dtime,dalpha)
        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesLambda4_y (dx,dy,dtime,dalpha)
        case (5)
          ! 5.)
          Dvalues(1) = fct_stokesLambda5_y (dx,dy,dtime,dalpha)
        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
      
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
