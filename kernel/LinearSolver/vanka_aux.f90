!##############################################################################
!# ****************************************************************************
!# <name> vanka_aux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of auxiliary routines for the VANKA solvers,
!# mainly highly-tuned routines for solving local systems.
!# </purpose>
!##############################################################################

module vanka_aux

use fsystem
use mprimitives

implicit none

contains

  ! ***************************************************************************

!<subroutine>

  pure subroutine vanka_aux_solve_BS2D_js_414(Du,Df,Da1,Da2,Db1,Db2,&
                                              Dd1,Dd2,Dc,Dm1,Dm2,Dm3,Dn)
  
!<description>
! Solves a local 2D Boussinesq system, Jacobi-Style version:
! 4 DOFs per velocity
! 1 DOF per pressure
! 4 DOFs per temperature
!
! / A1 0  B1 M1 \
! | 0  A2 B2 M2 | * u = f
! | D1 D2 C  M3 |
! \ 0  0  0  N  /
!
! Where:
! A1/A2/N is a 4x4 diagonal matrix
! B1/B2 is a 4x1 matrix
! D1/D2/M3 is a 1x4 matrix
! M1/M2 is a 4x4 matrix
! C is a 1x1 diagonal matrix
!
!</description>
  
!<input>
  ! The local matrices.
  real(DP), dimension(4), intent(IN) :: Da1,Da2,Dn
  real(DP), dimension(4,1), intent(IN) :: Db1,Db2
  real(DP), dimension(1,4), intent(IN) :: Dd1,Dd2,Dm3
  real(DP), dimension(4,4), intent(IN) :: Dm1,Dm2
  real(DP), dimension(1), intent(IN) :: Dc
  
  ! The local RHS vector.
  real(DP), dimension(13), intent(IN) :: Df
!</input>
  
!<output>
  ! The local solution vector.
  real(DP), dimension(13), intent(OUT) :: Du
!</output>
!</subroutine>

  ! local variables
  real(DP), dimension(4) :: Di1, Di2, Dt1, Dt2, Dg1, Dg2
  real(DP) :: dS
  
    ! Invert A1 and A2
    Di1(1) = 1.0_DP / Da1(1)
    Di1(2) = 1.0_DP / Da1(2)
    Di1(3) = 1.0_DP / Da1(3)
    Di1(4) = 1.0_DP / Da1(4)
    Di2(1) = 1.0_DP / Da2(1)
    Di2(2) = 1.0_DP / Da2(2)
    Di2(3) = 1.0_DP / Da2(3)
    Di2(4) = 1.0_DP / Da2(4)

    ! Calculate temperature
    ! t := N^-1 * f_t
    Du(10) = Df(10) / Dn(1)
    Du(11) = Df(11) / Dn(2)
    Du(12) = Df(12) / Dn(3)
    Du(13) = Df(13) / Dn(4)
    
    ! Calculate new RHS 
    ! g_u := f_u - M*t
    Dg1(1) = Df(1)-Dm1(1,1)*Du(10)-Dm1(1,2)*Du(11)-Dm1(1,3)*Du(12)-Dm1(1,4)*Du(13)
    Dg1(2) = Df(2)-Dm1(2,1)*Du(10)-Dm1(2,2)*Du(11)-Dm1(2,3)*Du(12)-Dm1(2,4)*Du(13)
    Dg1(3) = Df(3)-Dm1(3,1)*Du(10)-Dm1(3,2)*Du(11)-Dm1(3,3)*Du(12)-Dm1(3,4)*Du(13)
    Dg1(4) = Df(4)-Dm1(4,1)*Du(10)-Dm1(4,2)*Du(11)-Dm1(4,3)*Du(12)-Dm1(4,4)*Du(13)
    Dg2(1) = Df(5)-Dm2(1,1)*Du(10)-Dm2(1,2)*Du(11)-Dm2(1,3)*Du(12)-Dm2(1,4)*Du(13)
    Dg2(2) = Df(6)-Dm2(2,1)*Du(10)-Dm2(2,2)*Du(11)-Dm2(2,3)*Du(12)-Dm2(2,4)*Du(13)
    Dg2(3) = Df(7)-Dm2(3,1)*Du(10)-Dm2(3,2)*Du(11)-Dm2(3,3)*Du(12)-Dm2(3,4)*Du(13)
    Dg2(4) = Df(8)-Dm2(4,1)*Du(10)-Dm2(4,2)*Du(11)-Dm2(4,3)*Du(12)-Dm2(4,4)*Du(13)
    
    ! Precalculate D * A^-1
    Dt1(1) = Dd1(1,1)*Di1(1)
    Dt1(2) = Dd1(1,2)*Di1(2)
    Dt1(3) = Dd1(1,3)*Di1(3)
    Dt1(4) = Dd1(1,4)*Di1(4)
    Dt2(1) = Dd2(1,1)*Di2(1)
    Dt2(2) = Dd2(1,2)*Di2(2)
    Dt2(3) = Dd2(1,3)*Di2(3)
    Dt2(4) = Dd2(1,4)*Di2(4)

    ! Calculate Schur-Complement of A
    ! S := -C + D * A^-1 * B 
    dS = -Dc(1) &
       + Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1) &
       + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)
    
    ! Calculate pressure
    ! p := D * A^-1 * g_u + M3 * t - f_p
    Du(9) = (-Df(9) &
       + Dm3(1,1)*Du(10) + Dm3(1,2)*Du(11) &
       + Dm3(1,3)*Du(12) + Dm3(1,4)*Du(13) &
       + Dt1(1)*Dg1(1)+Dt1(2)*Dg1(2)+Dt1(3)*Dg1(3)+Dt1(4)*Dg1(4) &
       + Dt2(1)*Dg2(1)+Dt2(2)*Dg2(2)+Dt2(3)*Dg2(3)+Dt2(4)*Dg2(4)) / dS
    
    ! Calculate X- and Y-velocity
    ! u := A^-1 * (g_u - B * p)
    Du(1) = di1(1)*(Dg1(1) - Db1(1,1)*Du(9))
    Du(2) = di1(2)*(Dg1(2) - Db1(2,1)*Du(9))
    Du(3) = di1(3)*(Dg1(3) - Db1(3,1)*Du(9))
    Du(4) = di1(4)*(Dg1(4) - Db1(4,1)*Du(9))
    Du(5) = di2(1)*(Dg2(1) - Db2(1,1)*Du(9))
    Du(6) = di2(2)*(Dg2(2) - Db2(2,1)*Du(9))
    Du(7) = di2(3)*(Dg2(3) - Db2(3,1)*Du(9))
    Du(8) = di2(4)*(Dg2(4) - Db2(4,1)*Du(9))
    
    ! That's it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine vanka_aux_solve_BS2D_bd_414(Du,Df,Da1,Da2,Db1,Db2,&
                                              Dd1,Dd2,Dc,Dm1,Dm2,Dm3,Dn)
  
!<description>
! Solves a local 2D Boussinesq system, Block-Diagonal version:
! 4 DOFs per velocity
! 1 DOF per pressure
! 4 DOFs per temperature
!
! / A1 0  B1 M1 \
! | 0  A2 B2 M2 | * u = f
! | D1 D2 C  M3 |
! \ 0  0  0  N  /
!
! Where:
! A1/A2/N/M1/M2 is a 4x4 matrix
! B1/B2 is a 4x1 matrix
! D1/D2/M3 is a 1x4 matrix
! C is a 1x1 matrix
!
!</description>
  
!<input>
  real(DP), dimension(4,1), intent(IN) :: Db1,Db2
  real(DP), dimension(1,4), intent(IN) :: Dd1,Dd2,Dm3
  real(DP), dimension(4,4), intent(IN) :: Da1,Da2,Dm1,Dm2,Dn
  real(DP), dimension(13), intent(IN) :: Df
  real(DP), dimension(1,1), intent(IN) :: Dc
!</input>
  
!<output>
  real(DP), dimension(13), intent(OUT) :: Du
!</output>
!</subroutine>

  ! local variables
  real(DP), dimension(4,4) :: Di1, Di2, Dj
  real(DP), dimension(4) :: Dt1,Dt2,Dg1,Dg2
  real(DP) :: dS
  
    ! Invert A1, A2 and N
    call mprim_invert4x4MatrixDirectDble(Da1, Di1)
    call mprim_invert4x4MatrixDirectDble(Da2, Di2)
    call mprim_invert4x4MatrixDirectDble(Dn, Dj)

    ! Calculate temperature
    ! t := N^-1 * f_t
    Du(10) = Dj(1,1)*Df(10)+Dj(1,2)*Df(11)+Dj(1,3)*Df(12)+Dj(1,4)*Df(13)
    Du(11) = Dj(2,1)*Df(10)+Dj(2,2)*Df(11)+Dj(2,3)*Df(12)+Dj(2,4)*Df(13)
    Du(12) = Dj(3,1)*Df(10)+Dj(3,2)*Df(11)+Dj(3,3)*Df(12)+Dj(3,4)*Df(13)
    Du(13) = Dj(4,1)*Df(10)+Dj(4,2)*Df(11)+Dj(4,3)*Df(12)+Dj(4,4)*Df(13)
    
    ! Calculate new RHS
    ! g_u := f_u - M*t
    Dg1(1) = Df(1)-Dm1(1,1)*Du(10)-Dm1(1,2)*Du(11)-Dm1(1,3)*Du(12)-Dm1(1,4)*Du(13)
    Dg1(2) = Df(2)-Dm1(2,1)*Du(10)-Dm1(2,2)*Du(11)-Dm1(2,3)*Du(12)-Dm1(2,4)*Du(13)
    Dg1(3) = Df(3)-Dm1(3,1)*Du(10)-Dm1(3,2)*Du(11)-Dm1(3,3)*Du(12)-Dm1(3,4)*Du(13)
    Dg1(4) = Df(4)-Dm1(4,1)*Du(10)-Dm1(4,2)*Du(11)-Dm1(4,3)*Du(12)-Dm1(4,4)*Du(13)
    Dg2(1) = Df(5)-Dm2(1,1)*Du(10)-Dm2(1,2)*Du(11)-Dm2(1,3)*Du(12)-Dm2(1,4)*Du(13)
    Dg2(2) = Df(6)-Dm2(2,1)*Du(10)-Dm2(2,2)*Du(11)-Dm2(2,3)*Du(12)-Dm2(2,4)*Du(13)
    Dg2(3) = Df(7)-Dm2(3,1)*Du(10)-Dm2(3,2)*Du(11)-Dm2(3,3)*Du(12)-Dm2(3,4)*Du(13)
    Dg2(4) = Df(8)-Dm2(4,1)*Du(10)-Dm2(4,2)*Du(11)-Dm2(4,3)*Du(12)-Dm2(4,4)*Du(13)
    
    ! Precalculate D * A^-1
    Dt1(1) = Dd1(1,1)*Di1(1,1)+Dd1(1,2)*Di1(2,1)+Dd1(1,3)*Di1(3,1)+Dd1(1,4)*Di1(4,1)
    Dt1(2) = Dd1(1,1)*Di1(1,2)+Dd1(1,2)*Di1(2,2)+Dd1(1,3)*Di1(3,2)+Dd1(1,4)*Di1(4,2)
    Dt1(3) = Dd1(1,1)*Di1(1,3)+Dd1(1,2)*Di1(2,3)+Dd1(1,3)*Di1(3,3)+Dd1(1,4)*Di1(4,3)
    Dt1(4) = Dd1(1,1)*Di1(1,4)+Dd1(1,2)*Di1(2,4)+Dd1(1,3)*Di1(3,4)+Dd1(1,4)*Di1(4,4)
    Dt2(1) = Dd2(1,1)*Di2(1,1)+Dd2(1,2)*Di2(2,1)+Dd2(1,3)*Di2(3,1)+Dd2(1,4)*Di2(4,1)
    Dt2(2) = Dd2(1,1)*Di2(1,2)+Dd2(1,2)*Di2(2,2)+Dd2(1,3)*Di2(3,2)+Dd2(1,4)*Di2(4,2)
    Dt2(3) = Dd2(1,1)*Di2(1,3)+Dd2(1,2)*Di2(2,3)+Dd2(1,3)*Di2(3,3)+Dd2(1,4)*Di2(4,3)
    Dt2(4) = Dd2(1,1)*Di2(1,4)+Dd2(1,2)*Di2(2,4)+Dd2(1,3)*Di2(3,4)+Dd2(1,4)*Di2(4,4)
    
    ! Calculate Schur-Complement of A
    ! S := -C + D * A^-1 * B 
    dS = -Dc(1,1)+Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1)&
                 +Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)
    
    ! Calculate pressure
    ! p := D * A^-1 * g_u + M3 * t - f_p
    Du(9) = (-Df(9) &
          + Dt1(1)*Dg1(1)+Dt1(2)*Dg1(2)+Dt1(3)*Dg1(3)+Dt1(4)*Dg1(4) &
          + Dt2(1)*Dg2(1)+Dt2(2)*Dg2(2)+Dt2(3)*Dg2(3)+Dt2(4)*Dg2(4) &
          + Dm3(1,1)*Du(10) + Dm3(1,2)*Du(11) &
          + Dm3(1,3)*Du(12) + Dm3(1,4)*Du(13)) / dS

    ! Update RHS
    ! g_u := g_u - B * p
    Dg1(1) = Dg1(1) - Db1(1,1)*Du(9)
    Dg1(2) = Dg1(2) - Db1(2,1)*Du(9)
    Dg1(3) = Dg1(3) - Db1(3,1)*Du(9)
    Dg1(4) = Dg1(4) - Db1(4,1)*Du(9)
    Dg2(1) = Dg2(1) - Db2(1,1)*Du(9)
    Dg2(2) = Dg2(2) - Db2(2,1)*Du(9)
    Dg2(3) = Dg2(3) - Db2(3,1)*Du(9)
    Dg2(4) = Dg2(4) - Db2(4,1)*Du(9)
    
    ! Calculate X- and Y-velocity
    ! u := A^-1 * g_u
    Du(1) = Di1(1,1)*Dg1(1)+Di1(1,2)*Dg1(2)+Di1(1,3)*Dg1(3)+Di1(1,4)*Dg1(4)
    Du(2) = Di1(2,1)*Dg1(1)+Di1(2,2)*Dg1(2)+Di1(2,3)*Dg1(3)+Di1(2,4)*Dg1(4)
    Du(3) = Di1(3,1)*Dg1(1)+Di1(3,2)*Dg1(2)+Di1(3,3)*Dg1(3)+Di1(3,4)*Dg1(4)
    Du(4) = Di1(4,1)*Dg1(1)+Di1(4,2)*Dg1(2)+Di1(4,3)*Dg1(3)+Di1(4,4)*Dg1(4)
    Du(5) = Di2(1,1)*Dg2(1)+Di2(1,2)*Dg2(2)+Di2(1,3)*Dg2(3)+Di2(1,4)*Dg2(4)
    Du(6) = Di2(2,1)*Dg2(1)+Di2(2,2)*Dg2(2)+Di2(2,3)*Dg2(3)+Di2(2,4)*Dg2(4)
    Du(7) = Di2(3,1)*Dg2(1)+Di2(3,2)*Dg2(2)+Di2(3,3)*Dg2(3)+Di2(3,4)*Dg2(4)
    Du(8) = Di2(4,1)*Dg2(1)+Di2(4,2)*Dg2(2)+Di2(4,3)*Dg2(3)+Di2(4,4)*Dg2(4)
    
    ! That's it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine vanka_aux_solve_NS3D_js_61(Du,Df,Da1,Da2,Da3,Db1,Db2,Db3,&
                                             Dd1,Dd2,Dd3,Dc)
  
!<description>
! Solves a local 3D Navier-Stokes system, Jacobi-Style version:
! 6 DOFs per velocity
! 1 DOF per pressure
!
! / A1 0  0  B1 \
! | 0  A2 0  B2 | * u = f
! | 0  0  A3 B3 |
! \ D1 D2 D3 C  /
!
! Where:
! A1/A2/A3 is a 6x6 diagonal matrix
! B1/B2/B3 is a 6x1 matrix
! D1/D2/D3 is a 1x4 matrix
! C is a 1x1 diagonal matrix
!
!</description>
  
!<input>
  ! The local matrices.
  real(DP), dimension(6), intent(IN) :: Da1,Da2,Da3
  real(DP), dimension(6,1), intent(IN) :: Db1,Db2,Db3
  real(DP), dimension(1,6), intent(IN) :: Dd1,Dd2,Dd3
  real(DP), dimension(1), intent(IN) :: Dc
  
  ! The local RHS vector.
  real(DP), dimension(19), intent(IN) :: Df
!</input>
  
!<output>
  ! The local solution vector.
  real(DP), dimension(19), intent(OUT) :: Du
!</output>
!</subroutine>

  ! local variables
  real(DP), dimension(6) :: Di1, Di2, Di3, Dt1,Dt2,Dt3
  real(DP) :: dS
  
    ! Invert A1, A2 and A3
    Di1(1) = 1.0_DP / Da1(1)
    Di1(2) = 1.0_DP / Da1(2)
    Di1(3) = 1.0_DP / Da1(3)
    Di1(4) = 1.0_DP / Da1(4)
    Di1(5) = 1.0_DP / Da1(5)
    Di1(5) = 1.0_DP / Da1(6)
    Di2(1) = 1.0_DP / Da2(1)
    Di2(2) = 1.0_DP / Da2(2)
    Di2(3) = 1.0_DP / Da2(3)
    Di2(4) = 1.0_DP / Da2(4)
    Di2(5) = 1.0_DP / Da2(5)
    Di2(6) = 1.0_DP / Da2(6)
    Di3(1) = 1.0_DP / Da3(1)
    Di3(2) = 1.0_DP / Da3(2)
    Di3(3) = 1.0_DP / Da3(3)
    Di3(4) = 1.0_DP / Da3(4)
    Di3(5) = 1.0_DP / Da3(5)
    Di3(6) = 1.0_DP / Da3(6)
    
    ! Precalculate D * A^-1
    Dt1(1) = Dd1(1,1)*Di1(1)
    Dt1(2) = Dd1(1,2)*Di1(2)
    Dt1(3) = Dd1(1,3)*Di1(3)
    Dt1(4) = Dd1(1,4)*Di1(4)
    Dt1(5) = Dd1(1,5)*Di1(5)
    Dt1(6) = Dd1(1,6)*Di1(6)
    Dt2(1) = Dd2(1,1)*Di2(1)
    Dt2(2) = Dd2(1,2)*Di2(2)
    Dt2(3) = Dd2(1,3)*Di2(3)
    Dt2(4) = Dd2(1,4)*Di2(4)
    Dt2(5) = Dd2(1,5)*Di2(5)
    Dt2(6) = Dd2(1,6)*Di2(6)
    Dt3(1) = Dd3(1,1)*Di3(1)
    Dt3(2) = Dd3(1,2)*Di3(2)
    Dt3(3) = Dd3(1,3)*Di3(3)
    Dt3(4) = Dd3(1,4)*Di3(4)
    Dt3(5) = Dd3(1,5)*Di3(5)
    Dt3(6) = Dd3(1,6)*Di3(6)

    ! Calculate Schur-Complement of A
    ! S := -C + D * A^-1 * B 
    dS = -Dc(1) &
       + Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1) &
       + Dt1(4)*Db1(4,1)+Dt1(5)*Db1(5,1)+Dt1(6)*Db1(6,1) &
       + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1) &
       + Dt2(4)*Db2(4,1)+Dt2(5)*Db2(5,1)+Dt2(6)*Db2(6,1) &
       + Dt3(1)*Db3(1,1)+Dt3(2)*Db3(2,1)+Dt3(3)*Db3(3,1) &
       + Dt3(4)*Db3(4,1)+Dt3(5)*Db3(5,1)+Dt3(6)*Db3(6,1)
    
    ! Calculate pressure
    ! p := D * A^-1 * g_u - f_p
    Du(19) = (-Df(19) &
          + Dt1(1)*Df( 1)+Dt1(2)*Df( 2)+Dt1(3)*Df( 3) &
          + Dt1(4)*Df( 4)+Dt1(5)*Df( 5)+Dt1(6)*Df( 6) &
          + Dt2(1)*Df( 7)+Dt2(2)*Df( 8)+Dt2(3)*Df( 9) &
          + Dt2(4)*Df(10)+Dt2(5)*Df(11)+Dt2(6)*Df(12) &
          + Dt3(1)*Df(13)+Dt3(2)*Df(14)+Dt3(3)*Df(15) &
          + Dt3(4)*Df(16)+Dt3(5)*Df(17)+Dt3(6)*Df(18)) / dS
    
    ! Calculate X-, Y- and Z-velocity
    ! u := A^-1 * (f_u - B * p)
    Du( 1) = di1(1)*(Df( 1) - Db1(1,1)*Du(19))
    Du( 2) = di1(2)*(Df( 2) - Db1(2,1)*Du(19))
    Du( 3) = di1(3)*(Df( 3) - Db1(3,1)*Du(19))
    Du( 4) = di1(4)*(Df( 4) - Db1(4,1)*Du(19))
    Du( 5) = di1(5)*(Df( 5) - Db1(5,1)*Du(19))
    Du( 6) = di1(6)*(Df( 6) - Db1(6,1)*Du(19))
    Du( 7) = di2(1)*(Df( 7) - Db2(1,1)*Du(19))
    Du( 8) = di2(2)*(Df( 8) - Db2(2,1)*Du(19))
    Du( 9) = di2(3)*(Df( 9) - Db2(3,1)*Du(19))
    Du(10) = di2(4)*(Df(10) - Db2(4,1)*Du(19))
    Du(11) = di2(5)*(Df(11) - Db2(5,1)*Du(19))
    Du(12) = di2(6)*(Df(12) - Db2(6,1)*Du(19))
    Du(13) = di3(1)*(Df(13) - Db3(1,1)*Du(19))
    Du(14) = di3(2)*(Df(14) - Db3(2,1)*Du(19))
    Du(15) = di3(3)*(Df(15) - Db3(3,1)*Du(19))
    Du(16) = di3(4)*(Df(16) - Db3(4,1)*Du(19))
    Du(17) = di3(5)*(Df(17) - Db3(5,1)*Du(19))
    Du(18) = di3(6)*(Df(18) - Db3(6,1)*Du(19))
    
    ! That's it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine vanka_aux_solve_NS3D_bd_61(Du,Df,Da1,Da2,Da3,Db1,Db2,Db3,&
                                             Dd1,Dd2,Dd3,Dc)
  
!<description>
! Solves a local 3D Navier-Stokes system, block-diagonal version:
! 6 DOFs per velocity
! 1 DOF per pressure
!
! / A1 0  0  B1 \
! | 0  A2 0  B2 | * u = f
! | 0  0  A3 B3 |
! \ D1 D2 D3 C  /
!
! Where:
! A1/A2/A3 is a 6x6 matrix
! B1/B2/B3 is a 6x1 matrix
! D1/D2/D3 is a 1x6 matrix
! C is a 1x1 matrix
!
!</description>
  
!<input>
  real(DP), dimension(6,1), intent(IN) :: Db1,Db2,Db3
  real(DP), dimension(1,6), intent(IN) :: Dd1,Dd2,Dd3
  real(DP), dimension(6,6), intent(IN) :: Da1,Da2,Da3
  real(DP), dimension(19), intent(IN) :: Df
  real(DP), dimension(1,1), intent(IN) :: Dc
!</input>
  
!<output>
  real(DP), dimension(19), intent(OUT) :: Du
!</output>
!</subroutine>

  ! local variables
  real(DP), dimension(6,6) :: Di1, Di2, Di3
  real(DP), dimension(6) :: Dt1,Dt2,Dt3,Dg1,Dg2,Dg3
  real(DP) :: dS
  
    ! Invert A1, A2 and A3
    call mprim_invert6x6MatrixDirectDble(Da1, Di1)
    call mprim_invert6x6MatrixDirectDble(Da2, Di2)
    call mprim_invert6x6MatrixDirectDble(Da3, Di3)

    ! Precalculate D * A^-1
    Dt1(1) = Dd1(1,1)*Di1(1,1)+Dd1(1,2)*Di1(2,1)+Dd1(1,3)*Di1(3,1) &
           + Dd1(1,4)*Di1(4,1)+Dd1(1,5)*Di1(5,1)+Dd1(1,6)*Di1(6,1)
    Dt1(2) = Dd1(1,1)*Di1(1,2)+Dd1(1,2)*Di1(2,2)+Dd1(1,3)*Di1(3,2) &
           + Dd1(1,4)*Di1(4,2)+Dd1(1,5)*Di1(5,2)+Dd1(1,6)*Di1(6,2)
    Dt1(3) = Dd1(1,1)*Di1(1,3)+Dd1(1,2)*Di1(2,3)+Dd1(1,3)*Di1(3,3) &
           + Dd1(1,4)*Di1(4,3)+Dd1(1,5)*Di1(5,3)+Dd1(1,6)*Di1(6,3)
    Dt1(4) = Dd1(1,1)*Di1(1,4)+Dd1(1,2)*Di1(2,4)+Dd1(1,3)*Di1(3,4) &
           + Dd1(1,4)*Di1(4,4)+Dd1(1,5)*Di1(5,4)+Dd1(1,6)*Di1(6,4)
    Dt1(5) = Dd1(1,1)*Di1(1,5)+Dd1(1,2)*Di1(2,5)+Dd1(1,3)*Di1(3,5) &
           + Dd1(1,4)*Di1(4,5)+Dd1(1,5)*Di1(5,5)+Dd1(1,6)*Di1(6,5)
    Dt1(6) = Dd1(1,1)*Di1(1,6)+Dd1(1,2)*Di1(2,6)+Dd1(1,3)*Di1(3,6) &
           + Dd1(1,4)*Di1(4,6)+Dd1(1,5)*Di1(5,6)+Dd1(1,6)*Di1(6,6)
    Dt2(1) = Dd2(1,1)*Di2(1,1)+Dd2(1,2)*Di2(2,1)+Dd2(1,3)*Di2(3,1) &
           + Dd2(1,4)*Di2(4,1)+Dd2(1,5)*Di2(5,1)+Dd2(1,6)*Di2(6,1)
    Dt2(2) = Dd2(1,1)*Di2(1,2)+Dd2(1,2)*Di2(2,2)+Dd2(1,3)*Di2(3,2) &
           + Dd2(1,4)*Di2(4,2)+Dd2(1,5)*Di2(5,2)+Dd2(1,6)*Di2(6,2)
    Dt2(3) = Dd2(1,1)*Di2(1,3)+Dd2(1,2)*Di2(2,3)+Dd2(1,3)*Di2(3,3) &
           + Dd2(1,4)*Di2(4,3)+Dd2(1,5)*Di2(5,3)+Dd2(1,6)*Di2(6,3)
    Dt2(4) = Dd2(1,1)*Di2(1,4)+Dd2(1,2)*Di2(2,4)+Dd2(1,3)*Di2(3,4) &
           + Dd2(1,4)*Di2(4,4)+Dd2(1,5)*Di2(5,4)+Dd2(1,6)*Di2(6,4)
    Dt2(5) = Dd2(1,1)*Di2(1,5)+Dd2(1,2)*Di2(2,5)+Dd2(1,3)*Di2(3,5) &
           + Dd2(1,4)*Di2(4,5)+Dd2(1,5)*Di2(5,5)+Dd2(1,6)*Di2(6,5)
    Dt2(6) = Dd2(1,1)*Di2(1,6)+Dd2(1,2)*Di2(2,6)+Dd2(1,3)*Di2(3,6) &
           + Dd2(1,4)*Di2(4,6)+Dd2(1,5)*Di2(5,6)+Dd2(1,6)*Di2(6,6)
    Dt3(1) = Dd3(1,1)*Di3(1,1)+Dd3(1,2)*Di3(2,1)+Dd3(1,3)*Di3(3,1) &
           + Dd3(1,4)*Di3(4,1)+Dd3(1,5)*Di3(5,1)+Dd3(1,6)*Di3(6,1)
    Dt3(2) = Dd3(1,1)*Di3(1,2)+Dd3(1,2)*Di3(2,2)+Dd3(1,3)*Di3(3,2) &
           + Dd3(1,4)*Di3(4,2)+Dd3(1,5)*Di3(5,2)+Dd3(1,6)*Di3(6,2)
    Dt3(3) = Dd3(1,1)*Di3(1,3)+Dd3(1,2)*Di3(2,3)+Dd3(1,3)*Di3(3,3) &
           + Dd3(1,4)*Di3(4,3)+Dd3(1,5)*Di3(5,3)+Dd3(1,6)*Di3(6,3)
    Dt3(4) = Dd3(1,1)*Di3(1,4)+Dd3(1,2)*Di3(2,4)+Dd3(1,3)*Di3(3,4) &
           + Dd3(1,4)*Di3(4,4)+Dd3(1,5)*Di3(5,4)+Dd3(1,6)*Di3(6,4)
    Dt3(5) = Dd3(1,1)*Di3(1,5)+Dd3(1,2)*Di3(2,5)+Dd3(1,3)*Di3(3,5) &
           + Dd3(1,4)*Di3(4,5)+Dd3(1,5)*Di3(5,5)+Dd3(1,6)*Di3(6,5)
    Dt3(6) = Dd3(1,1)*Di3(1,6)+Dd3(1,2)*Di3(2,6)+Dd3(1,3)*Di3(3,6) &
           + Dd3(1,4)*Di3(4,6)+Dd3(1,5)*Di3(5,6)+Dd3(1,6)*Di3(6,6)
    
    ! Calculate Schur-Complement of A
    ! S := -C + D * A^-1 * B 
    dS = -Dc(1,1) &
       + Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1) &
       + Dt1(4)*Db1(4,1)+Dt1(5)*Db1(5,1)+Dt1(6)*Db1(6,1) &
       + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1) &
       + Dt2(4)*Db2(4,1)+Dt2(5)*Db2(5,1)+Dt2(6)*Db2(6,1) &
       + Dt3(1)*Db3(1,1)+Dt3(2)*Db3(2,1)+Dt3(3)*Db3(3,1) &
       + Dt3(4)*Db3(4,1)+Dt3(5)*Db3(5,1)+Dt3(6)*Db3(6,1)
    
    ! Calculate pressure
    ! p := D * A^-1 * g_u - f_p
    Du(19) = (-Df(19) &
          + Dt1(1)*Df( 1)+Dt1(2)*Df( 2)+Dt1(3)*Df( 3) &
          + Dt1(4)*Df( 4)+Dt1(5)*Df( 5)+Dt1(6)*Df( 6) &
          + Dt2(1)*Df( 7)+Dt2(2)*Df( 8)+Dt2(3)*Df( 9) &
          + Dt2(4)*Df(10)+Dt2(5)*Df(11)+Dt2(6)*Df(12) &
          + Dt3(1)*Df(13)+Dt3(2)*Df(14)+Dt3(3)*Df(15) &
          + Dt3(4)*Df(16)+Dt3(5)*Df(17)+Dt3(6)*Df(18)) / dS

    ! Update RHS
    ! g_u := f_u - B * p
    Dg1(1) = Df( 1) - Db1(1,1)*Du(19)
    Dg1(2) = Df( 2) - Db1(2,1)*Du(19)
    Dg1(3) = Df( 3) - Db1(3,1)*Du(19)
    Dg1(4) = Df( 4) - Db1(4,1)*Du(19)
    Dg1(5) = Df( 5) - Db1(5,1)*Du(19)
    Dg1(6) = Df( 6) - Db1(6,1)*Du(19)
    Dg2(1) = Df( 7) - Db2(1,1)*Du(19)
    Dg2(2) = Df( 8) - Db2(2,1)*Du(19)
    Dg2(3) = Df( 9) - Db2(3,1)*Du(19)
    Dg2(4) = Df(10) - Db2(4,1)*Du(19)
    Dg2(5) = Df(11) - Db2(5,1)*Du(19)
    Dg2(6) = Df(12) - Db2(6,1)*Du(19)
    Dg3(1) = Df(13) - Db3(1,1)*Du(19)
    Dg3(2) = Df(14) - Db3(2,1)*Du(19)
    Dg3(3) = Df(15) - Db3(3,1)*Du(19)
    Dg3(4) = Df(16) - Db3(4,1)*Du(19)
    Dg3(5) = Df(17) - Db3(5,1)*Du(19)
    Dg3(6) = Df(18) - Db3(6,1)*Du(19)
    
    ! Calculate X-, Y- and Z-velocity
    ! u := A^-1 * g_u
    Du( 1) = Di1(1,1)*Dg1(1)+Di1(1,2)*Dg1(2)+Di1(1,3)*Dg1(3)&
           + Di1(1,4)*Dg1(4)+Di1(1,5)*Dg1(5)+Di1(1,6)*Dg1(6)
    Du( 2) = Di1(2,1)*Dg1(1)+Di1(2,2)*Dg1(2)+Di1(2,3)*Dg1(3)&
           + Di1(2,4)*Dg1(4)+Di1(2,5)*Dg1(5)+Di1(2,6)*Dg1(6)
    Du( 3) = Di1(3,1)*Dg1(1)+Di1(3,2)*Dg1(2)+Di1(3,3)*Dg1(3)&
           + Di1(3,4)*Dg1(4)+Di1(3,5)*Dg1(5)+Di1(3,6)*Dg1(6)
    Du( 4) = Di1(4,1)*Dg1(1)+Di1(4,2)*Dg1(2)+Di1(4,3)*Dg1(3)&
           + Di1(4,4)*Dg1(4)+Di1(4,5)*Dg1(5)+Di1(4,6)*Dg1(6)
    Du( 5) = Di1(5,1)*Dg1(1)+Di1(5,2)*Dg1(2)+Di1(5,3)*Dg1(3)&
           + Di1(5,4)*Dg1(4)+Di1(5,5)*Dg1(5)+Di1(5,6)*Dg1(6)
    Du( 6) = Di1(6,1)*Dg1(1)+Di1(6,2)*Dg1(2)+Di1(6,3)*Dg1(3)&
           + Di1(6,4)*Dg1(4)+Di1(6,5)*Dg1(5)+Di1(6,6)*Dg1(6)
    Du( 7) = Di2(1,1)*Dg2(1)+Di2(1,2)*Dg2(2)+Di2(1,3)*Dg2(3)&
           + Di2(1,4)*Dg2(4)+Di2(1,5)*Dg2(5)+Di2(1,6)*Dg2(6)
    Du( 8) = Di2(2,1)*Dg2(1)+Di2(2,2)*Dg2(2)+Di2(2,3)*Dg2(3)&
           + Di2(2,4)*Dg2(4)+Di2(2,5)*Dg2(5)+Di2(2,6)*Dg2(6)
    Du( 9) = Di2(3,1)*Dg2(1)+Di2(3,2)*Dg2(2)+Di2(3,3)*Dg2(3)&
           + Di2(3,4)*Dg2(4)+Di2(3,5)*Dg2(5)+Di2(3,6)*Dg2(6)
    Du(10) = Di2(4,1)*Dg2(1)+Di2(4,2)*Dg2(2)+Di2(4,3)*Dg2(3)&
           + Di2(4,4)*Dg2(4)+Di2(4,5)*Dg2(5)+Di2(4,6)*Dg2(6)
    Du(11) = Di2(5,1)*Dg2(1)+Di2(5,2)*Dg2(2)+Di2(5,3)*Dg2(3)&
           + Di2(5,4)*Dg2(4)+Di2(5,5)*Dg2(5)+Di2(5,6)*Dg2(6)
    Du(12) = Di2(6,1)*Dg2(1)+Di2(6,2)*Dg2(2)+Di2(6,3)*Dg2(3)&
           + Di2(6,4)*Dg2(4)+Di2(6,5)*Dg2(5)+Di2(6,6)*Dg2(6)
    Du(13) = Di3(1,1)*Dg3(1)+Di3(1,2)*Dg3(2)+Di3(1,3)*Dg3(3)&
           + Di3(1,4)*Dg3(4)+Di3(1,5)*Dg3(5)+Di3(1,6)*Dg3(6)
    Du(14) = Di3(2,1)*Dg3(1)+Di3(2,2)*Dg3(2)+Di3(2,3)*Dg3(3)&
           + Di3(2,4)*Dg3(4)+Di3(2,5)*Dg3(5)+Di3(2,6)*Dg3(6)
    Du(15) = Di3(3,1)*Dg3(1)+Di3(3,2)*Dg3(2)+Di3(3,3)*Dg3(3)&
           + Di3(3,4)*Dg3(4)+Di3(3,5)*Dg3(5)+Di3(3,6)*Dg3(6)
    Du(16) = Di3(4,1)*Dg3(1)+Di3(4,2)*Dg3(2)+Di3(4,3)*Dg3(3)&
           + Di3(4,4)*Dg3(4)+Di3(4,5)*Dg3(5)+Di3(4,6)*Dg3(6)
    Du(17) = Di3(5,1)*Dg3(1)+Di3(5,2)*Dg3(2)+Di3(5,3)*Dg3(3)&
           + Di3(5,4)*Dg3(4)+Di3(5,5)*Dg3(5)+Di3(5,6)*Dg3(6)
    Du(18) = Di3(6,1)*Dg3(1)+Di3(6,2)*Dg3(2)+Di3(6,3)*Dg3(3)&
           + Di3(6,4)*Dg3(4)+Di3(6,5)*Dg3(5)+Di3(6,6)*Dg3(6)
    
    ! That's it
  
  end subroutine

end module