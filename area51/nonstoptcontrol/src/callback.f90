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
  use derivatives
  
  use scalarpde
  use linearformevaluation
  use collection
  use domainintegration
  use paramlist
  
  use physics
  
  implicit none

contains

  ! ***************************************************************************

  subroutine cb_getPhysics(rparlist,rphysics)
  
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
        "dpar", rphysics%dpar)
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
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "dtimemin", rphysics%dtimemin)
    call parlst_getvalue_double (rparlist, "PHYSICS", &
        "dtimemax", rphysics%dtimemax)
  
  end subroutine

  ! ***************************************************************************
  ! Reference functions, Heat equation
  ! ***************************************************************************
  ! ***************************************************************************
  ! Heat equation, function set 1.
  !   y = t^2   x1 
  ! ***************************************************************************

  elemental real(DP) function fct_heatY1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY1 = dtime**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda1 = dalpha * (- 2.0_DP*dtime * dx) + 2*dalpha*dx
  end function
    
  elemental real(DP) function fct_heatF1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF1 = 2.0_DP*dtime*dx - 2.0_DP*dx*(dtime-1.0_DP)
  end function
  
  elemental real(DP) function fct_heatZ1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatZ1 = -2.0_DP*dalpha*dx + dtime**2*dx
  end function

  ! ***************************************************************************
  ! Heat equation, function set 2
  !   y = t^2 (1-t)^2   x1 
  !   l = a t (1-t) x1
  ! No factor 4 due to inconsistent terminal condition.
  ! ***************************************************************************
  
  elemental real(DP) function fct_heatY2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY2 = dtime**2 * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda2 = dalpha * dtime * (1.0_DP-dtime) * dx 
  end function

  elemental real(DP) function fct_heatF2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF2 = 3.0_DP*dtime*dx-7.0_DP*dtime**2*dx+4*dtime**3*dx
  end function

  elemental real(DP) function fct_heatZ2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatZ2 = (1.0_DP-dtime)*dx-dtime*dx+dtime**2*(1.0_DP-dtime)**2*dx 
  end function

  ! ***************************************************************************
  ! Heat equation, function set 3
  !   y = t^2 (1-t)^2   x1 
  !   l = a t (1-t)^2   x1
  ! Gives factor 4 with CN due to proper terminal condition.
  ! ***************************************************************************

  elemental real(DP) function fct_heatY3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY3 = dtime**2 * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda3 = dalpha * dtime * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatF3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF3 = &
      3.0_DP*dtime*dx-8.0_DP*dtime**2*dx+5*dtime**3*dx
  end function

  elemental real(DP) function fct_heatZ3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatZ3 = &
         (1.0_DP-dtime)**2*dx-2.0_DP*dtime*(1.0_DP-dtime)*dx+&
          dtime**2*(1.0_DP-dtime)**2*dx 
  end function
          
  ! ***************************************************************************
  ! Heat equation, function set 4
  !   y = t^2 (1-t)^2   x1 
  ! ***************************************************************************
    
  elemental real(DP) function fct_heatY4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY4 = dtime**2 * (1.0_DP-dtime)**2 * dx 
  end function

  elemental real(DP) function fct_heatLambda4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda4 = dalpha * (-2.0_DP*dtime*(1.0_DP-dtime)**2*dx + 2.0_DP*dtime**2*(1.0_DP-dtime)*dx)
  end function
    
  elemental real(DP) function fct_heatF4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF4 = 0.0_DP
  end function

  elemental real(DP) function fct_heatZ4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
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

  elemental real(DP) function fct_heatY5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY5 = dtime * dx 
  end function

  elemental real(DP) function fct_heatLambda5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda5 = (1.0_DP-dtime)*dx
  end function
    
  elemental real(DP) function fct_heatF5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF5 = -(dtime-dalpha-1.0_DP)*dx/dalpha
  end function
  
  elemental real(DP) function fct_heatZ5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatZ5 = (dtime-1.0_DP)*dx
  end function

  ! ***************************************************************************
  ! Heat equation, function set 6.
  !   y = t^2   x1 
  !   l = (1-t)^2 x1
  ! ***************************************************************************

  elemental real(DP) function fct_heatY6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY6 = dtime**2 * dx
  end function

  elemental real(DP) function fct_heatLambda6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda6 = (1.0_DP-dtime)**2*dx
  end function
    
  elemental real(DP) function fct_heatF6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF6 = dx*(2.0_DP*dtime*dalpha+1.0_DP-2.0_DP*dtime+dtime**2)/dalpha
  end function
  
  elemental real(DP) function fct_heatZ6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatZ6 = dtime**2*dx-2.0_DP*(1.0_DP-dtime)*dx
  end function

  ! ***************************************************************************
  ! Heat equation, function set 7.
  !   wa = exp(da*PI^2*t) sin(PI*x1) sin(PI*x2)
  !    y = -1(2+da) PI^2 wa(t,x1,x2)
  !    l = wa(t,x1,x2) - wa(1.0,x1,x2)
  ! ***************************************************************************

  elemental real(DP) function fct_eig7 (da,dx,dy,dtime,dalpha)
  real(DP), intent(in) :: da,dx,dy,dtime,dalpha
    fct_eig7 = exp(da*SYS_PI**2*dtime) * sin(SYS_PI*dx) * sin(SYS_PI*dy)
  end function

  elemental real(DP) function fct_heatY7 (da,dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: da,dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY7 = -SYS_PI**2/(2.0_DP+da) * fct_eig7 (da,dx,dy,dtime,dalpha)
  end function

  elemental real(DP) function fct_heatLambda7 (da,dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: da,dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda7 = fct_eig7 (da,dx,dy,dtime,dalpha) - fct_eig7 (da,dx,dy,dtimemax,dalpha)
  end function
    
  elemental real(DP) function fct_heatF7 (da,dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: da,dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF7 = - ( SYS_PI**4 * dalpha * fct_eig7 (da,dx,dy,dtime,dalpha) - &
                     fct_eig7 (da,dx,dy,dtime,dalpha) + &
                     fct_eig7 (da,dx,dy,dtimemax,dalpha) ) / dalpha
    !fct_heatF7 = -SYS_PI**4 * fct_eig7 (da,dx,dy,1.0_DP,dalpha)
  end function
  
  elemental real(DP) function fct_heatZ7 (da,dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: da,dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatZ7 = SYS_PI**2 / (2.0_DP+da) * &
        ( da**2 * fct_eig7 (da,dx,dy,dtime,dalpha) - &
          5.0_DP * fct_eig7 (da,dx,dy,dtime,dalpha) + &
          4.0_DP * fct_eig7 (da,dx,dy,dtimemax,dalpha) + &
          2.0_DP * da * fct_eig7 (da,dx,dy,dtimemax,dalpha) )
    !fct_heatZ7 = (da**2-2.0_DP)/(1.0_DP+da) * SYS_PI**2 * fct_eig7 (da,dx,dy,dtime,dalpha) + &
    !    SYS_PI**2*fct_eig7 (da,dx,dy,1.0_DP,dalpha)
  end function

  ! ***************************************************************************
  ! Heat equation, function set 8.
  !   w = sin(Pi*x1)*sin(Pi*x2)
  !
  !   y = sin(PI*t/2) * w
  !   l = (sin(PI*t/2)-1) * w
  !
  !   f = (1/2)*w*(cos((1/2)*Pi*t)*Pi*alpha+4*alpha*Pi^2*sin((1/2)*Pi*t)+2*sin((1/2)*Pi*t)-2)/alpha
  !   z = -(1/2)*w*(-cos((1/2)*Pi*t)*Pi+4*Pi^2*sin((1/2)*Pi*t)-4*Pi^2-2*sin((1/2)*Pi*t))
  ! ***************************************************************************

  elemental real(DP) function fct_eig8 (dx,dy,dtime,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dalpha
    fct_eig8 = sin(SYS_PI*dx) * sin(SYS_PI*dy)
  end function

  elemental real(DP) function fct_heatY8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatY8 = sin(0.5_DP*SYS_PI*dtime) * fct_eig8 (dx,dy,dtime,dalpha)
  end function

  elemental real(DP) function fct_heatLambda8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatLambda8 = fct_eig8 (dx,dy,dtime,dalpha) * (sin(0.5_DP*SYS_PI*dtime) - 1.0_DP)
  end function
    
  elemental real(DP) function fct_heatF8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatF8 = 0.5_DP/dalpha * fct_eig8 (dx,dy,dtime,dalpha) * &
      ( SYS_PI*dalpha*cos(0.5_DP*SYS_PI*dtime) + &
        4.0_DP*dalpha*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime) + &
        2.0_DP*sin(0.5_DP*SYS_PI*dtime) - 2.0_DP )
  end function
  
  elemental real(DP) function fct_heatZ8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_heatZ8 = -0.5_DP * fct_eig8 (dx,dy,dtime,dalpha) * &
      ( -SYS_PI*cos(0.5_DP*SYS_PI*dtime) + 4.0_DP*SYS_PI**2*sin(0.5_DP*SYS_PI*dtime) & 
        - 4.0_DP*SYS_PI**2 - 2.0_DP*sin(0.5_DP*SYS_PI*dtime) )
  end function

  ! ***************************************************************************
  ! Reference functions, Stokes equation
  ! ***************************************************************************
  ! ***************************************************************************
  ! Stokes, function set 1.
  !   y = t^2  ( x1 , -x2 , 0)
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY1_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY1_x = dtime**2 * dx
  end function

  elemental real(DP) function fct_stokesY1_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY1_y = - dtime**2 * dy
  end function

  elemental real(DP) function fct_stokesP1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP1 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda1_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda1_x = dalpha * (-2.0_DP*dtime*dx)
  end function

  elemental real(DP) function fct_stokesLambda1_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda1_y = dalpha * (2.0_DP*dtime*dy)
  end function

  elemental real(DP) function fct_stokesXi1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi1 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ1_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ1_x = -2.0_DP*dx + dtime**2*dx
  end function

  elemental real(DP) function fct_stokesZ1_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ1_y = 2.0_DP*dy - dtime**2*dy
  end function

  ! ***************************************************************************
  ! Stokes, function set 2.
  !   y = t^2 (1-t)^3  ( x2(1-x2) , 0 , 0 )
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY2_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY2_x = dy * (-1 + dy) * dtime ** 2 * (-1 + dtime) ** 3
  end function

  elemental real(DP) function fct_stokesY2_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY2_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP2 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda2_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda2_x = &
        - dalpha * dtime * (-1 + dtime) ** 2 * (2 * dtime - 2 * dtime ** 2 &
        + 2 * dy - 5 * dy * dtime - 2 * dy ** 2 + 5 * dy ** 2 * dtime)
  end function

  elemental real(DP) function fct_stokesLambda2_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda2_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi2 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ2_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ2_x = &
          dy * dtime ** 2 - 3 * dy * dtime ** 3 + 3 * dy * dtime ** 4 &
        - dy * dtime ** 5 - dy ** 2 * dtime ** 2 + 3 * dy ** 2 * dtime ** 3 &
        - 3 * dy ** 2 * dtime ** 4 + dy ** 2 * dtime ** 5 - 2 * dalpha * dy &
        - 0.18D2 * dalpha * dy ** 2 * dtime + 2 * dalpha * dy ** 2 &
        + 0.18D2 * dalpha * dy * dtime - 0.36D2 * dalpha * dy * dtime ** 2 &
        + 0.20D2 * dalpha * dy * dtime ** 3 + 0.36D2 * dalpha * dy ** 2 * dtime ** 2 &
        - 0.20D2 * dalpha * dy ** 2 * dtime ** 3
  end function

  elemental real(DP) function fct_stokesZ2_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ2_y = 0.0_DP
  end function

  ! ***************************************************************************
  ! Stokes, function set 3.
  !   y = t^2 (x1^2 x2 , -x1 x2^2 , 0)
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY3_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY3_x = dtime**2 * dx**2 * dy
  end function

  elemental real(DP) function fct_stokesY3_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY3_y = - dtime**2 * dx**2 * dy
  end function

  elemental real(DP) function fct_stokesP3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP3 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda3_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda3_x = -dalpha*(-2*dy*dtime**2*(1-dtime)**2 &
        +2*dx**2*dy*dtime*(1-dtime)**2 &
        -2*dx**2*dy*dtime**2*(1-dtime) )
  end function

  elemental real(DP) function fct_stokesLambda3_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda3_y = -dalpha*(-2*dx*dtime**2*(1-dtime)**2 &
        +2*dx*dy**2*dtime*(1-dtime)**2 &
        -2*dx*dy**2*dtime**2*(1-dtime) )
  end function

  elemental real(DP) function fct_stokesXi3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi3 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ3_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ3_x = (dx**2*dy*dtime**2*(-1+dtime)**2) &
        - dalpha*(-4*dy*dtime+12*dy*dtime**2 &
                  -8*dy*dtime**3+2*dx**2*dy &
                  -12*dx**2*dy*dtime+12*dx**2*dy*dtime**2)
  end function

  elemental real(DP) function fct_stokesZ3_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ3_y = (dx*dy**2*dtime**2*(-1+dtime)**2) &
          - dalpha*(-4*dx*dtime+12*dx*dtime**2 &
                    -8*dx*dtime**3+2*dx*dy**2 &
                    -12*dx*dy**2*dtime+12*dx*dy**2*dtime**2)
  end function

  ! ***************************************************************************
  ! Stokes, function set 4.
  !   y = t^2 (1-t)^2  (x1^3 x2^2 , -x1^2 x2^3 , 0)
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY4_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY4_x = dx ** 3 * dy ** 2 * dtime ** 2 * (-1.0_DP + dtime) ** 2
  end function

  elemental real(DP) function fct_stokesY4_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY4_y = -dx ** 2 * dy ** 3 * dtime ** 2 * (-1.0_DP + dtime) ** 2
  end function

  elemental real(DP) function fct_stokesP4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP4 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda4_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda4_x = -2 * dalpha * dx * dtime * (-1 + dtime) * &
              (3 * dy ** 2 * dtime - 3 * dy ** 2 * dtime ** 2 + dx ** 2 * dtime &
              - dx ** 2 * dtime ** 2 - dx ** 2 * dy ** 2 + 2 * dx ** 2 * dy ** 2 * dtime)
  end function

  elemental real(DP) function fct_stokesLambda4_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda4_y = 2 * dalpha * dy * dtime * (-1 + dtime) * &
              (dy ** 2 * dtime - dy ** 2 * dtime ** 2 + 3.0_DP * dx ** 2 * dtime &
              - 3.0_DP * dx ** 2 * dtime ** 2 - dx ** 2 * dy ** 2 + 2 * dx ** 2 * dy ** 2 * dtime)
  end function

  elemental real(DP) function fct_stokesXi4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi4 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ4_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ4_x = &
        -2*dalpha*dx**3*dy**2+12*dalpha*dx**3*dy**2*dtime-12*dalpha*dx**3*dy**2*dtime**2+dx**3*dy**2*dtime**2 &
        -2*dx**3*dy**2*dtime**3+dx**3*dy**2*dtime**4+24*dalpha*dx*dtime**2-48*dalpha*dx*dtime**3+24*dalpha*dx*dtime**4
  end function

  elemental real(DP) function fct_stokesZ4_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ4_y = &
        2*dalpha*dy**3*dx**2-12*dalpha*dy**3*dx**2*dtime+12*dalpha*dy**3*dx**2*dtime**2-dx**2*dy**3*dtime**2 &
        +2*dx**2*dy**3*dtime**3-dx**2*dy**3*dtime**4-24*dalpha*dy*dtime**2+48*dalpha*dy*dtime**3-24*dalpha*dy*dtime**4
  end function
  
  ! ***************************************************************************
  ! Stokes, function set 5.
  !   y = t  ( x2(1-x2) , 0 , 0 )
  ! No coupling between primal and dual!
  ! EXAMPLE DOES NOT WORK, NOT DIVERGENCE FREE WHEN BEING TIME DEPENDENT!
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY5_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY5_x = dtime * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesY5_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY5_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP5 = dtime*(1.0_DP-dx) !0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda5_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda5_x = 0.0_DP !-dalpha*( &
!                2.0_DP*dtime**2*(1.0_DP-dtime)**2 &
!              + 2.0_DP*dy*(1.0_DP-dy)*dtime*(1.0_DP-dtime)**2 &
!              - 2.0_DP*dy*(1.0_DP-dy)*dtime**2*(1.0_DP-dtime))
  end function

  elemental real(DP) function fct_stokesLambda5_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda5_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi5 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ5_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ5_x = 0.0_DP !- dy*(-1.0_DP+dy)*dtime**2*(-1.0_DP+dtime)**2 &
!            - dalpha*(4.0_DP*dtime-12.0_DP*dtime**2+8.0_DP*dtime**3+2.0_DP*dy &
!                      -12.0_DP*dy*dtime+12.0_DP*dy*dtime**2 &
!                      -2.0_DP*dy**2+12.0_DP*dy**2*dtime &
!                      -12.0_DP*dy**2*dtime**2)
  end function

  elemental real(DP) function fct_stokesZ5_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ5_y = 0.0_DP
  end function

  ! ***************************************************************************
  ! Stokes, function set 6.
  !   y      =     t ( x2(1-x2) , 0 , 0 )
  !   lambda = (1-t) ( x2(1-x2) , 0 , 0 )
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY6_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY6_x = dtime * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesY6_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY6_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP6 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda6_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda6_x = (1.0_DP-dtime) * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesLambda6_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda6_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi6 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ6_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ6_x = -0.2D1 + 0.2D1 * dtime - dy + dy ** 2 + dy * dtime &
        - dy ** 2 * dtime
  end function

  elemental real(DP) function fct_stokesZ6_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ6_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesF6_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF6_x = (0.2D1 * dtime * dalpha + dy * dalpha - dy ** 2 * dalpha + dy &
        - dy * dtime - dy ** 2 + dy ** 2 * dtime) / dalpha
  end function

  elemental real(DP) function fct_stokesF6_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF6_y = 0.0_DP
  end function

  ! ***************************************************************************
  ! Stokes, function set 7.
  !   y      =     t^2 ( x2(1-x2) , 0 , 0 )
  !   lambda = (1-t)^2 ( x2(1-x2) , 0 , 0 )
  ! ***************************************************************************

  elemental real(DP) function fct_stokesY7_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY7_x = dtime**2 * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesY7_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesY7_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesP7 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesP7 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesLambda7_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda7_x = (1.0_DP-dtime)**2 * dy*(1.0_DP-dy)
  end function

  elemental real(DP) function fct_stokesLambda7_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesLambda7_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesXi7 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesXi7 = 0.0_DP
  end function

  elemental real(DP) function fct_stokesZ7_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ7_x = -0.2D1 + 0.4D1 * dtime - 0.2D1 * dtime ** 2 &
        - 0.2D1 * dy + 0.2D1 * dy * dtime + 0.2D1 * dy ** 2 &
        - 0.2D1 * dy ** 2 * dtime + dy * dtime ** 2 - dy ** 2 * dtime ** 2
  end function

  elemental real(DP) function fct_stokesZ7_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesZ7_y = 0.0_DP
  end function

  elemental real(DP) function fct_stokesF7_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
    fct_stokesF7_x = -(-0.2D1 * dtime ** 2 * dalpha &
        - 0.2D1 * dy * dtime * dalpha + 0.2D1 * dy ** 2 * dtime * dalpha &
        - dy + 0.2D1 * dy * dtime - dy * dtime ** 2 + dy ** 2 &
        - 0.2D1 * dy ** 2 * dtime + dy ** 2 * dtime ** 2) / dalpha
  end function

  elemental real(DP) function fct_stokesF7_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
  real(DP), intent(in) :: dx,dy,dtime,dtimemin,dtimemax,dalpha
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

  ! ***************************************************************************

!<subroutine>

  subroutine ferrFunction (cderivative, rdiscretisation, &
                                  nelements, npointsPerElement, Dpoints, &
                                  IdofsTest, rdomainIntSubset, &
                                  Dvalues, rcollection)
  
  use fsystem
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
!<description>
  ! Callback function to calculate the error of a FE function to the
  ! corresponding analytical function.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in) :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
  
!</subroutine>
    integer :: icomponent,cequation,creferenceProblem,ierrType
    real(DP) :: dtime,dalpha,dpar,dtimemin,dtimemax

    ! Get the settings of the equation
    icomponent = rcollection%IquickAccess (1)
    cequation = rcollection%IquickAccess (2)
    creferenceProblem = rcollection%IquickAccess (3)
    dtime = rcollection%DquickAccess (1)
    dalpha = rcollection%DquickAccess (2)
    dpar = rcollection%DquickAccess (3)
    dtimemin = rcollection%DquickAccess (4)
    dtimemax = rcollection%DquickAccess (5)

    if (cderivative .eq. DER_FUNC) then
      ! Select the equation and calculate the error.
      select case (cequation)
      case (0,2)
        ! -------------------------------------------------------------
        ! Heat equation
        ! -------------------------------------------------------------
        select case (icomponent)
        case (1)
          !Dvalues(:,:) = dtime*(Dpoints(1,:,:)*Dpoints(2,:,:))
          !Dvalues(:,:) = 0.0_DP
          
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY1 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY1 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
          
          case (2)
            ! 2.) -> BC in spacetimebc.f90 beachten!
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY2 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY2 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case (3)
            ! 3.) -> BC in spacetimebc.f90 beachten!
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY3 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY3 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case (4)
            ! 4.) -> BC in spacetimebc.f90 beachten!
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY4 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY4 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case (5)
            ! 5.) -> BC in spacetimebc.f90 beachten!
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY5 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY5 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case (6)
            ! 6.) -> BC in spacetimebc.f90 beachten!
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY6 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY6 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case (7)
            ! 7.) -> BC in spacetimebc.f90 beachten!
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY7 (dpar,Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY7 (dpar,Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case (8)
            ! 8.) -> BC in spacetimebc.f90 beachten!
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatY8 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatY8 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
          
        case (2)
          !Dvalues(:,:) = (1.0_DP-dtime)*(Dpoints(1,:,:)*Dpoints(2,:,:))
          !Dvalues(:,:) = 0.0_DP
          !Dvalues(:,:) = - (2.0_DP*dtime**2 * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) 
          ! +  2.0_DP*dtime**2.0_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )
          
          ! For 1D: See BC in spacetimebc.f90!!!
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda1 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda1 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
            
          case (2)
            ! 2.) 
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda2 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda2 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if

          case (3)
            ! 3.) 
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda3 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda3 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
            
          case (4)
            ! 4.) 
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda4 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda4 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
            
          case (5)
            ! 5.) 
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda5 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda5 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
            
          case (6)
            ! 6.) 
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda6 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda6 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
            
          case (7)
            ! 7.) 
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda7 (dpar,Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda7 (dpar,Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
            
          case (8)
            ! 8.) 
            if (ubound(Dpoints,1) .eq. NDIM1D) then
              Dvalues(:,:) = fct_heatLambda8 (Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha)
            else
              Dvalues(:,:) = fct_heatLambda8 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
            end if
            
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
            Dvalues(:,:) = fct_stokesY1_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
          
          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesY2_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
          
          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesY3_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesY4_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesY5_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesY6_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesY7_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (8)
            ! 8.)
            Dvalues(:,:) = fct_stokesY8_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (9)
            ! 9.)
            Dvalues(:,:) = fct_stokesY9_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
          
        case (2)
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesY1_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesY2_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesY3_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesY4_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesY5_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesY6_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesY7_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (8)
            ! 8.)
            Dvalues(:,:) = fct_stokesY8_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (9)
            ! 9.)
            Dvalues(:,:) = fct_stokesY9_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
        
        case (3)
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesP1 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesP2 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesP3 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesP4 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesP5 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesP6 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesP7 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (8)
            ! 8.)
            Dvalues(:,:) = fct_stokesP8 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (9)
            ! 9.)
            Dvalues(:,:) = fct_stokesP9 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

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
            Dvalues(:,:) = fct_stokesLambda1_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesLambda2_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesLambda3_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesLambda4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesLambda5_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesLambda6_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesLambda7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (8)
            ! 8.)
            Dvalues(:,:) = fct_stokesLambda8_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (9)
            ! 9.)
            Dvalues(:,:) = fct_stokesLambda9_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
          
        case (5)
        
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesLambda1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
          
          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesLambda2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
          
          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesLambda3_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesLambda4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesLambda5_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesLambda6_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesLambda7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (8)
            ! 8.)
            Dvalues(:,:) = fct_stokesLambda8_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (9)
            ! 9.)
            Dvalues(:,:) = fct_stokesLambda9_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
        
        case (6)
        
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesXi1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesXi2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesXi3(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesXi4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesXi5(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesXi6(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesXi7(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (8)
            ! 8.)
            Dvalues(:,:) = fct_stokesXi8(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

          case (9)
            ! 9.)
            Dvalues(:,:) = fct_stokesXi9(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)

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

    else
    
      call output_line ("only function values supported.")
      call sys_halt()
      
    end if

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
  real(DP) :: dcouplePrimalToDual,dcoupleDualToPrimal
  real(DP) :: dtimeMin,dtimeMax

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  
    integer :: icomponent,creferenceProblem,ithetaschemetype
    real(DP) :: dtime,dalpha,dtstep,dgamma,dtheta,dpar
    
    icomponent = rcollection%IquickAccess(1)
    creferenceProblem = rcollection%IquickAccess(2)
    ithetaschemetype = rcollection%IquickAccess(3)
    dtime = rcollection%DquickAccess(1)
    dtstep = rcollection%DquickAccess(2)
    dalpha = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    dtheta = rcollection%DquickAccess(5)
    dpar = rcollection%DquickAccess(6)
    dcouplePrimalToDual = rcollection%DquickAccess(7)
    dcoupleDualToPrimal = rcollection%DquickAccess(8)
    dtimeMin = rcollection%DquickAccess(9)
    dtimeMax = rcollection%DquickAccess(10)
    
    Dcoefficients(1,:,:) = 0.0_DP
    
    ! Primal or dual?
    select case (icomponent)
    case (1)
      ! Primal: f

      !Dcoefficients(1,:,:) = Dpoints(1,:,:)*Dpoints(2,:,:) &
      !                     + ((1.0_DP-dtime)/dalpha)*Dpoints(1,:,:)*Dpoints(2,:,:)
      !Dcoefficients(1,:,:) = Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
      !                     + 2*dtime*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))) 
      !                     + ((1.0_DP-dtime)/dalpha)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
      
      !Dcoefficients(1,:,:) = 2*dtime*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
      !                     + 2*dtime**2*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)))        
      !                     + ((1.0_DP-dtime)**2/dalpha)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
      
      select case (creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF1(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatF1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if
      
      case (2)
        ! 2.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF2(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatF2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (3)
        ! 3.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF3(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatF3(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (4)
        ! 4.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF4(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatF4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (5)
        ! 5.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF5(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatF5(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (6)
        ! 6.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF6(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatF6(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (7)
        ! 7.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF7(dpar,Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = &
            fct_heatF7(dpar,Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (8)
        ! 8.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatF8(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = &
            fct_heatF8(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select
                            
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
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ1(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if
      
      case (2)  
        ! 2.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ2(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (3)  
        ! 3.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ3(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ3(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if
              
      case (4)
        ! 4.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ4(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (5)
        ! 5.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ5(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ5(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (6)
        ! 6.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ6(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ6(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (7)
        ! 7.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ7(dpar,Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ7(dpar,Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (8)
        ! 8.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            - (fct_heatZ8(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            - (fct_heatZ8(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select
      
!      if (ithetaschemetype .eq. 1) then
!        ! Dual solution lives in the interval midpoints.
!        ! Dual RHS lives in the endpoints of the interval!
!        ! dtime = 1.0 is ok!
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (dtheta + dtheta*dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP-dtheta)
!        end if
!      else
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP + dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = 0.0_DP
!        end if
!      end if
    case default
      ! Should not happen
      call sys_halt()
    end select      

  end subroutine

  ! ***************************************************************************

  subroutine sol_heatEquation (rdiscretisation, rform, &
                nelements, npointsPerElement, Dpoints, &
                IdofsTest, rdomainIntSubset, &
                Dcoefficients, rcollection)
  
  ! Coefficient function for the solution of the heat equation
  ! if being used in the RHS.
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  type(t_linearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset    
  type(t_collection), intent(inout), optional :: rcollection
  real(DP) :: dcouplePrimalToDual,dcoupleDualToPrimal
  real(DP) :: dtimeMin,dtimeMax

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  
    integer :: icomponent,creferenceProblem,ithetaschemetype
    real(DP) :: dtime,dalpha,dtstep,dgamma,dtheta,dpar
    
    icomponent = rcollection%IquickAccess(1)
    creferenceProblem = rcollection%IquickAccess(2)
    ithetaschemetype = rcollection%IquickAccess(3)
    dtime = rcollection%DquickAccess(1)
    dtstep = rcollection%DquickAccess(2)
    dalpha = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    dtheta = rcollection%DquickAccess(5)
    dpar = rcollection%DquickAccess(6)
    dcouplePrimalToDual = rcollection%DquickAccess(7)
    dcoupleDualToPrimal = rcollection%DquickAccess(8)
    dtimeMin = rcollection%DquickAccess(9)
    dtimeMax = rcollection%DquickAccess(10)
    
    Dcoefficients(1,:,:) = 0.0_DP
    
    ! Primal or dual?
    select case (icomponent)
    case (1)
      ! Primal: y

      !Dcoefficients(1,:,:) = Dpoints(1,:,:)*Dpoints(2,:,:) &
      !                     + ((1.0_DP-dtime)/dalpha)*Dpoints(1,:,:)*Dpoints(2,:,:)
      !Dcoefficients(1,:,:) = Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
      !                     + 2*dtime*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))) 
      !                     + ((1.0_DP-dtime)/dalpha)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
      
      !Dcoefficients(1,:,:) = 2*dtime*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) &
      !                     + 2*dtime**2*(Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)))        
      !                     + ((1.0_DP-dtime)**2/dalpha)*Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
      
      select case (creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY1(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatY1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if
      
      case (2)
        ! 2.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY2(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatY2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (3)
        ! 3.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY3(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatY3(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (4)
        ! 4.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY4(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatY4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (5)
        ! 5.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY5(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatY5(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (6)
        ! 6.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY6(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = fct_heatY6(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (7)
        ! 7.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY7(dpar,Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = &
            fct_heatY7(dpar,Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case (8)
        ! 8.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = fct_heatY8(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) 
        else
          Dcoefficients(1,:,:) = &
            fct_heatY8(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) 
        end if

      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select
                            
    case (2)
      ! Dual: lambda
      
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
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda1(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda1(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if
      
      case (2)  
        ! 2.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda2(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda2(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (3)  
        ! 3.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda3(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda3(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if
              
      case (4)
        ! 4.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda4(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda4(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (5)
        ! 5.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda5(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda5(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (6)
        ! 6.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda6(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda6(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (7)
        ! 7.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda7(dpar,Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda7(dpar,Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case (8)
        ! 8.)
        if (ubound(Dpoints,1) .eq. NDIM1D) then
          Dcoefficients(1,:,:) = &
            (fct_heatLambda8(Dpoints(1,:,:),0.0_DP,dtime,dtimeMin,dtimeMax,dalpha) )
        else
          Dcoefficients(1,:,:) = &
            (fct_heatLambda8(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
        end if

      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select
      
!      if (ithetaschemetype .eq. 1) then
!        ! Dual solution lives in the interval midpoints.
!        ! Dual RHS lives in the endpoints of the interval!
!        ! dtime = 1.0 is ok!
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (dtheta + dtheta*dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP-dtheta)
!        end if
!      else
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP + dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = 0.0_DP
!        end if
!      end if
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
  real(DP) :: dcouplePrimalToDual,dcoupleDualToPrimal

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  
    integer :: icomponent,creferenceProblem,ithetaschemetype
    real(DP) :: dtime,dalpha,dtstep,dtheta,dgamma,dtimeMin,dtimeMax
    
    icomponent = rcollection%IquickAccess(1)
    creferenceProblem = rcollection%IquickAccess(2)
    ithetaschemetype = rcollection%IquickAccess(3)
    dtime = rcollection%DquickAccess(1)
    dtstep = rcollection%DquickAccess(2)
    dalpha = rcollection%DquickAccess(3)
    dgamma = rcollection%DquickAccess(4)
    dtheta = rcollection%DquickAccess(5)
    dcouplePrimalToDual = rcollection%DquickAccess(7)
    dcoupleDualToPrimal = rcollection%DquickAccess(8)
    dtimeMin = rcollection%DquickAccess(9)
    dtimeMax = rcollection%DquickAccess(10)
    
    Dcoefficients(1,:,:) = 0.0_DP
    
    ! Primal or dual?
    select case (icomponent)
    case (1)
      ! Primal: f1
      select case (creferenceProblem)
      case (6)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_StokesLambda6_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF6_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      case (7)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_StokesLambda7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      case (8)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_StokesLambda8_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF8_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      case (9)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_stokesLambda9_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF9_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      end select

    case (2)
      ! Primal: f2
      select case (creferenceProblem)
      case (6)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_StokesLambda6_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF6_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      case (7)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_StokesLambda7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      case (8)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_StokesLambda8_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF8_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      case (9)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              -(1.0_DP-dcoupleDualToPrimal)/dalpha * &
                (fct_stokesLambda9_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
              + fct_stokesF9_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
      end select

    case (3)
      ! Primal: f_p
      select case (creferenceProblem)
      case (8)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              fct_stokesF8_p(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
        
      case (9)
        ! No RHS in initial condition since the solution should
        ! be =0 there.
        if (dtime .gt. 0.0_DP) then
          Dcoefficients(1,:,:) =  &
              fct_stokesF9_p(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)
        end if
        
      case default
        Dcoefficients(1,:,:) = 0.0_DP
      end select
      
    case (4)
      ! Dual: -z1
      
      select case (creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        Dcoefficients(1,:,:) =  &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY1_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ1_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (2)
        ! 2.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY2_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ2_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (3)
        ! 3.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY3_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ3_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      
      case (4)
        ! 4.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )

      case (5)
        ! 5.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY5_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ5_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      
      case (6)
        ! 6.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY6_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ6_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      
      case (7)
        ! 7.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      
      case (8)
        ! 8.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY8_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ8_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )

      case (9)
        ! 9.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY9_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ9_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )

      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select

!      if (ithetaschemetype .eq. 1) then
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (dtheta + dtheta*dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP-dtheta)
!        end if
!      else
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP + 1.0_DP*dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = 0.0_DP
!        end if
!      end if

    case (5)
      ! Dual: -z2
      
      select case (creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        Dcoefficients(1,:,:) =  &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ1_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (2)
        ! 2.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ2_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (3)
        ! 3.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY3_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ3_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (4)
        ! 4.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ4_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (5)
        ! 5.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY5_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ5_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (6)
        ! 6.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY6_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ6_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      case (7)
        ! 7.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ7_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )

      case (8)
        ! 8.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY8_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ8_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )

      case (9)
        ! 9.)
        Dcoefficients(1,:,:) = &
            (1.0_DP-dcouplePrimalToDual) * &
                (fct_stokesY9_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha)) &
            - (fct_stokesZ9_y(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )

      case default
        call output_line ("Problem not supported.")
        call sys_halt()
      end select

!      if (ithetaschemetype .eq. 1) then
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (dtheta + dtheta*dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP-dtheta)
!        end if
!      else
!        if (dtime .eq. 1.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = Dcoefficients(1,:,:) * (1.0_DP + dgamma/dtstep)
!        end if
!
!        if (dtime .eq. 0.0_DP) then
!          !Dcoefficients(1,:,:) = (-2.0_DP*Dpoints(1,:,:)/dtstep)
!          Dcoefficients(1,:,:) = 0.0_DP
!        end if
!      end if

    case (6)
      ! Dual: z_p
      select case (creferenceProblem)
      case (8)
        ! 8.)
        Dcoefficients(1,:,:) = &
            (fct_stokesZ8_p(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )

      case (9)
        ! 9.)
        Dcoefficients(1,:,:) = &
            (fct_stokesZ9_p(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dtimeMin,dtimeMax,dalpha) )
      
      case default
        Dcoefficients(1,:,:) = 0.0_DP
      end select

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
    real(DP) :: dtime,dalpha,dpar,dtimeMin,dtimeMax
    integer :: icomponent,cequation,creferenceProblem
    
    icomponent = rcollection%IquickAccess(1)
    cequation = rcollection%IquickAccess(2)
    creferenceProblem = rcollection%IquickAccess(3)
    dtime = rcollection%DquickAccess(1)
    dalpha = rcollection%DquickAccess(2)
    dpar = rcollection%DquickAccess(4)
    dtimeMin = rcollection%DquickAccess(5)
    dtimeMax = rcollection%DquickAccess(6)
    
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
          Dvalues(1) = fct_heatY1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (2)
          ! 2.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (3)
          ! 3.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (4)
          ! 4.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (5)
          ! 5.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (6)
          ! 6.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (7)
          ! 6.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY7 (dpar,dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (8)
          ! 8.) -> BC in spacetimebc.f90 beachten!
          Dvalues(1) = fct_heatY8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

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
          Dvalues(1) = fct_heatLambda1 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (2)
          ! 2.) 
          Dvalues(1) = fct_heatLambda2 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (3)
          ! 3.) 
          Dvalues(1) = fct_heatLambda3 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (4)
          ! 4.) 
          Dvalues(1) = fct_heatLambda4 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (5)
          ! 5.) 
          Dvalues(1) = fct_heatLambda5 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (6)
          ! 6.) 
          Dvalues(1) = fct_heatLambda6 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (7)
          ! 6.) 
          Dvalues(1) = fct_heatLambda7 (dpar,dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (8)
          ! 8.) 
          Dvalues(1) = fct_heatLambda8 (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
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
          Dvalues(1) = fct_stokesY1_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
    
        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesY2_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
     
        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesY3_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesY4_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (5)
          ! 5.)
          Dvalues(1) = fct_stokesY5_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (6)
          ! 6.)
          Dvalues(1) = fct_stokesY6_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (7)
          ! 7.)
          Dvalues(1) = fct_stokesY7_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (8)
          ! 8.)
          Dvalues(1) = fct_stokesY8_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (9)
          ! 9.)
          Dvalues(1) = fct_stokesY9_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
        
      case (2)
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_stokesY1_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesY2_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesY3_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesY4_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (5)
          ! 5.)
          Dvalues(1) = fct_stokesY5_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (6)
          ! 6.)
          Dvalues(1) = fct_stokesY6_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (7)
          ! 7.)
          Dvalues(1) = fct_stokesY7_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (8)
          ! 8.)
          Dvalues(1) = fct_stokesY8_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (9)
          ! 9.)
          Dvalues(1) = fct_stokesY9_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

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
          Dvalues(1) = fct_stokesLambda1_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesLambda2_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
          
        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesLambda3_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesLambda4_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (5)
          ! 5.)
          Dvalues(1) = fct_stokesLambda5_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (6)
          ! 6.)
          Dvalues(1) = fct_stokesLambda6_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (7)
          ! 7.)
          Dvalues(1) = fct_stokesLambda7_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (8)
          ! 8.)
          Dvalues(1) = fct_stokesLambda8_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case (9)
          ! 9.)
          Dvalues(1) = fct_stokesLambda9_x (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)

        case default
          call output_line ("Problem not supported.")
          call sys_halt()
        end select
        
      case (5)
      
        select case (creferenceProblem)
        case (0)
        case (1)
          ! 1.)
          Dvalues(1) = fct_stokesLambda1_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (2)
          ! 2.)
          Dvalues(1) = fct_stokesLambda2_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (3)
          ! 3.)
          Dvalues(1) = fct_stokesLambda3_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (4)
          ! 4.)
          Dvalues(1) = fct_stokesLambda4_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (5)
          ! 5.)
          Dvalues(1) = fct_stokesLambda5_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (6)
          ! 6.)
          Dvalues(1) = fct_stokesLambda6_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (7)
          ! 7.)
          Dvalues(1) = fct_stokesLambda7_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (8)
          ! 8.)
          Dvalues(1) = fct_stokesLambda8_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
        case (9)
          ! 9.)
          Dvalues(1) = fct_stokesLambda9_y (dx,dy,dtime,dtimeMin,dtimeMax,dalpha)
        
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
