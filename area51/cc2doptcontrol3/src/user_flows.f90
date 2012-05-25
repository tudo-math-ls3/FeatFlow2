!##############################################################################
!# ****************************************************************************
!# <name> user_flows </name>
!# ****************************************************************************
!#
!# <purpose>
!# Defines a couple of user-defined, analytically given flows and functions.
!# </purpose>
!##############################################################################

module user_flows

use fsystem

contains

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
    fct_stokesZ4_x = &
        -2*dalpha*dx**3*dy**2+12*dalpha*dx**3*dy**2*dtime-12*dalpha*dx**3*dy**2*dtime**2+dx**3*dy**2*dtime**2 &
        -2*dx**3*dy**2*dtime**3+dx**3*dy**2*dtime**4+24*dalpha*dx*dtime**2-48*dalpha*dx*dtime**3+24*dalpha*dx*dtime**4
!    fct_stokesZ4_x = 6 * dx * dy ** 2 * dtime ** 2 &
!        - 0.12D2 * dx * dy ** 2 * dtime ** 3 &
!        + 6 * dx * dy ** 2 * dtime ** 4 + 2 * dx ** 3 * dtime ** 2 &
!        - 4 * dx ** 3 * dtime ** 3 + 2 * dx ** 3 * dtime ** 4 &
!        + dx ** 3 * dy ** 2 * dtime ** 2 - 2 * dx ** 3 * dy ** 2 * dtime ** 3 &
!        + dx ** 3 * dy ** 2 * dtime ** 4 + 0.12D2 * dalpha * dx * dy ** 2 * dtime &
!        - 0.36D2 * dalpha * dx * dy ** 2 * dtime ** 2 &
!        + 0.24D2 * dalpha * dx * dy ** 2 * dtime ** 3 &
!        + 4 * dalpha * dx ** 3 * dtime &
!        - 0.12D2 * dalpha * dx ** 3 * dtime ** 2 &
!        + 8 * dalpha * dx ** 3 * dtime ** 3 &
!        - 2 * dalpha * dx ** 3 * dy ** 2 &
!        + 0.12D2 * dalpha * dx ** 3 * dy ** 2 * dtime &
!        - 0.12D2 * dalpha * dx ** 3 * dy ** 2 * dtime ** 2
  end function

  elemental real(DP) function fct_stokesZ4_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ4_y = &
        2*dalpha*dy**3*dx**2-12*dalpha*dy**3*dx**2*dtime+12*dalpha*dy**3*dx**2*dtime**2-dx**2*dy**3*dtime**2 &
        +2*dx**2*dy**3*dtime**3-dx**2*dy**3*dtime**4-24*dalpha*dy*dtime**2+48*dalpha*dy*dtime**3-24*dalpha*dy*dtime**4
!    fct_stokesZ4_y = -2 * dy ** 3 * dtime ** 2 + 4 * dy ** 3 * dtime ** 3 &
!        - 2 * dy ** 3 * dtime ** 4 - 6 * dx ** 2 * dy * dtime ** 2 &
!        + 0.12D2 * dx ** 2 * dy * dtime ** 3 - 6 * dx ** 2 * dy * dtime ** 4 &
!        - dx ** 2 * dy ** 3 * dtime ** 2 + 2 * dx ** 2 * dy ** 3 * dtime ** 3 &
!        - dx ** 2 * dy ** 3 * dtime ** 4 - 4 * dalpha * dy ** 3 * dtime &
!        + 0.12D2 * dalpha * dy ** 3 * dtime ** 2 &
!        - 8 * dalpha * dy ** 3 * dtime ** 3 &
!        - 0.12D2 * dalpha * dx ** 2 * dy * dtime &
!        + 0.36D2 * dalpha * dx ** 2 * dy * dtime ** 2 &
!        - 0.24D2 * dalpha * dx ** 2 * dy * dtime ** 3 &
!        + 2 * dalpha * dx ** 2 * dy ** 3 &
!        - 0.12D2 * dalpha * dx ** 2 * dy ** 3 * dtime &
!        + 0.12D2 * dalpha * dx ** 2 * dy ** 3 * dtime ** 2
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
  ! Stokes, function set 6.
  !   y      =     t ( x2(1-x2) , 0 , 0 )
  !   lambda = (1-t) ( x2(1-x2) , 0 , 0 )
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY6_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY6_x = dtime * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesY6_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY6_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP6 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesP6 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda6_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda6_x = (1.0_DP-dtime) * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesLambda6_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda6_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi6 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesXi6 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ6_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ6_x = -0.2D1 + 0.2D1 * dtime - dy + dy ** 2 + dy * dtime &
        - dy ** 2 * dtime
  end function

  elemental real(DP) function fct_stokesZ6_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ6_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesF6_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesF6_x = (0.2D1 * dtime * dalpha + dy * dalpha - dy ** 2 * dalpha + dy &
        - dy * dtime - dy ** 2 + dy ** 2 * dtime) / dalpha
  end function

  elemental real(DP) function fct_stokesF6_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesF6_y = 0.0_DP
  end function

  ! ***************************************************************************
  ! Stokes, function set 7.
  !   y      =     t^2 ( x2(1-x2) , 0 , 0 )
  !   lambda = (1-t)^2 ( x2(1-x2) , 0 , 0 )
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY7_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY7_x = dtime**2 * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesY7_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesY7_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP7 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesP7 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda7_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda7_x = (1.0_DP-dtime)**2 * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesLambda7_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesLambda7_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi7 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesXi7 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ7_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ7_x = -0.2D1 + 0.4D1 * dtime - 0.2D1 * dtime ** 2 &
        - 0.2D1 * dy + 0.2D1 * dy * dtime + 0.2D1 * dy ** 2 &
        - 0.2D1 * dy ** 2 * dtime + dy * dtime ** 2 - dy ** 2 * dtime ** 2
  end function

  elemental real(DP) function fct_stokesZ7_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesZ7_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesF7_x (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesF7_x = -(-0.2D1 * dtime ** 2 * dalpha &
        - 0.2D1 * dy * dtime * dalpha + 0.2D1 * dy ** 2 * dtime * dalpha &
        - dy + 0.2D1 * dy * dtime - dy * dtime ** 2 + dy ** 2 &
        - 0.2D1 * dy ** 2 * dtime + dy ** 2 * dtime ** 2) / dalpha
  end function

  elemental real(DP) function fct_stokesF7_y (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_stokesF7_y = 0.0_DP
  end function

  ! ***************************************************************************
  ! Stokes, function set 8
  !   w      = sin (Pi/2 x1) sin (Pi x2)
  !   y      = sin (Pi t ) w
  !   lambda = sin (Pi t ) w
  ! ***************************************************************************

  elemental real(DP) function fct_eigSt8(dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_eigSt8 = sin(0.5_DP*SYS_PI*dx) * sin(SYS_PI*dy)
  end function

  elemental real(DP) function fct_stokesY8_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY8_x = fct_eigSt8(dx,dy,dtime,dalpha)*sin(0.5_DP*SYS_PI*dtime)
  end function

  elemental real(DP) function fct_stokesY8_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY8_y = fct_eigSt8(dx,dy,dtime,dalpha)*sin(0.5_DP*SYS_PI*dtime)
  end function

  elemental real(DP) function fct_stokesP8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP8 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda8_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda8_x = fct_eigSt8(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1.0_DP)
  end function

  elemental real(DP) function fct_stokesLambda8_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda8_y = fct_eigSt8(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1.0_DP)
  end function

  elemental real(DP) function fct_stokesXi8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi8 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ8_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ8_x = -0.25_DP*fct_eigSt8(dx,dy,dtime,dalpha)*(-4.0_DP*sin(0.5_DP*SYS_PI*dtime)+&
        5.0_DP*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)-5.0_DP*SYS_PI**2-2*cos(0.5_DP*SYS_PI*dtime)*SYS_PI)
  end function

  elemental real(DP) function fct_stokesZ8_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ8_y = -0.25_DP*fct_eigSt8(dx,dy,dtime,dalpha)*(-4.0_DP*sin(0.5_DP*SYS_PI*dtime)+&
        5.0_DP*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)-5.0_DP*SYS_PI**2-2*cos(0.5_DP*SYS_PI*dtime)*SYS_PI)
  end function

  elemental real(DP) function fct_stokesZ8_p (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ8_p = -0.5_DP*SYS_PI*(sin(0.5_DP*SYS_PI*dtime)*cos(0.5_DP*SYS_PI*dx)*sin(SYS_PI*dy)+&
        2.0_DP*sin(0.5_DP*SYS_PI*dtime)*sin(0.5_DP*SYS_PI*dx)*cos(SYS_PI*dy)-&
        cos(0.5_DP*SYS_PI*dx)*sin(SYS_PI*dy)-2*sin(0.5_DP*SYS_PI*dx)*cos(SYS_PI*dy))
  end function

  elemental real(DP) function fct_stokesF8_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF8_x = 1.25_DP*fct_eigSt8(dx,dy,dtime,dalpha)*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)+&
        0.5_DP*fct_eigSt8(dx,dy,dtime,dalpha)*cos(0.5_DP*SYS_PI*dtime)*SYS_PI+&
        fct_eigSt8(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1)/dalpha
  end function

  elemental real(DP) function fct_stokesF8_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF8_y = 1.25_DP*fct_eigSt8(dx,dy,dtime,dalpha)*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)+&
        0.5_DP*fct_eigSt8(dx,dy,dtime,dalpha)*cos(0.5_DP*SYS_PI*dtime)*SYS_PI+&
        fct_eigSt8(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1)/dalpha
  end function

  elemental real(DP) function fct_stokesF8_p (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF8_p = -0.5_DP*SYS_PI*sin(0.5_DP*SYS_PI*dtime)*(cos(0.5_DP*SYS_PI*dx)*sin(SYS_PI*dy)+&
        2*sin(0.5_DP*SYS_PI*dx)*cos(SYS_PI*dy))
  end function


  ! ***************************************************************************
  ! Stokes, function set 9
  !   w      = sin (Pi/2 x1) sin (Pi x2)
  !   w_0    = sin (Pi x1) sin (Pi x2)
  !   y      = sin (Pi t ) ( w , w )
  !   lambda = sin (Pi t ) ( w , w )
  !   p      = sin (Pi t) w_0
  !   xi     = sin (Pi t) w_0
  ! ***************************************************************************

  elemental real(DP) function fct_eigSt9(dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_eigSt9 = sin(0.5_DP*SYS_PI*dx) * sin(SYS_PI*dy)
  end function

  elemental real(DP) function fct_eigSt9_1(dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    ! No int-mean = 0 since there is a Neumann edge that forces p=0 on the edge.
    fct_eigSt9_1 = sin(SYS_PI*dx) * sin(SYS_PI*dy)
  end function

  elemental real(DP) function fct_stokesY9_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY9_x = fct_eigSt9(dx,dy,dtime,dalpha)*sin(0.5_DP*SYS_PI*dtime)
  end function

  elemental real(DP) function fct_stokesY9_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY9_y = fct_eigSt9(dx,dy,dtime,dalpha)*sin(0.5_DP*SYS_PI*dtime)
  end function

  elemental real(DP) function fct_stokesP9 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP9 = sin(SYS_PI*dtime)*fct_eigSt9_1(dx,dy,dtime,dalpha)
  end function

  elemental real(DP) function fct_stokesLambda9_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda9_x = fct_eigSt9(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1.0_DP)
  end function

  elemental real(DP) function fct_stokesLambda9_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda9_y = fct_eigSt9(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1.0_DP)
  end function

  elemental real(DP) function fct_stokesXi9 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi9 = sin(SYS_PI*dtime)*fct_eigSt9_1(dx,dy,dtime,dalpha)
  end function

  elemental real(DP) function fct_stokesZ9_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ9_x = -0.25_DP*fct_eigSt9(dx,dy,dtime,dalpha)*(-4.0_DP*sin(0.5_DP*SYS_PI*dtime)+&
        5.0_DP*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)-5.0_DP*SYS_PI**2-2*cos(0.5_DP*SYS_PI*dtime)*SYS_PI) - &
        sin(SYS_PI*dtime)*SYS_PI*cos(SYS_PI*dx) * sin(SYS_PI*dy)
  end function

  elemental real(DP) function fct_stokesZ9_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ9_y = -0.25_DP*fct_eigSt9(dx,dy,dtime,dalpha)*(-4.0_DP*sin(0.5_DP*SYS_PI*dtime)+&
        5.0_DP*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)-5.0_DP*SYS_PI**2-2*cos(0.5_DP*SYS_PI*dtime)*SYS_PI) - &
        sin(SYS_PI*dtime)*SYS_PI*sin(SYS_PI*dx) * cos(SYS_PI*dy)
  end function

  elemental real(DP) function fct_stokesZ9_p (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ9_p = -0.5_DP*SYS_PI*(sin(0.5_DP*SYS_PI*dtime)*cos(0.5_DP*SYS_PI*dx)*sin(SYS_PI*dy)+&
        2.0_DP*sin(0.5_DP*SYS_PI*dtime)*sin(0.5_DP*SYS_PI*dx)*cos(SYS_PI*dy)-&
        cos(0.5_DP*SYS_PI*dx)*sin(SYS_PI*dy)-2*sin(0.5_DP*SYS_PI*dx)*cos(SYS_PI*dy))
  end function

  elemental real(DP) function fct_stokesF9_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF9_x = 1.25_DP*fct_eigSt9(dx,dy,dtime,dalpha)*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)+&
        0.5_DP*fct_eigSt9(dx,dy,dtime,dalpha)*cos(0.5_DP*SYS_PI*dtime)*SYS_PI+&
        fct_eigSt9(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1)/dalpha + &
        sin(SYS_PI*dtime)*SYS_PI*cos(SYS_PI*dx) * sin(SYS_PI*dy)
  end function

  elemental real(DP) function fct_stokesF9_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF9_y = 1.25_DP*fct_eigSt9(dx,dy,dtime,dalpha)*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime)+&
        0.5_DP*fct_eigSt9(dx,dy,dtime,dalpha)*cos(0.5_DP*SYS_PI*dtime)*SYS_PI+&
        fct_eigSt9(dx,dy,dtime,dalpha)*(sin(0.5_DP*SYS_PI*dtime)-1)/dalpha + &
        sin(SYS_PI*dtime)*SYS_PI*sin(SYS_PI*dx) * cos(SYS_PI*dy)
  end function

  elemental real(DP) function fct_stokesF9_p (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF9_p = -0.5_DP*SYS_PI*sin(0.5_DP*SYS_PI*dtime)*(cos(0.5_DP*SYS_PI*dx)*sin(SYS_PI*dy)+&
        2*sin(0.5_DP*SYS_PI*dx)*cos(SYS_PI*dy))
  end function

end module