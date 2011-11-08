!##############################################################################
!# ****************************************************************************
!# <name> optc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising a
!# stationary and nonstationary Navier-Stokes optimal control problem
!#
!#  $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alpha/2||u||^2 $$
!#
!#  $$- \nu Delta(y) + y*\Nabla(y) + \Nabla p = f $$
!#  $$ \Nabla \cdot y = 0$$
!#  $$- \nu Delta(\lambda) - y*\Nabla(\lambda) + (\Nabla y)^t\lambda + \Nabla \xi = y-z $$
!#  $$ \Nabla \cdot \lambda = 0$$
!#
!#
!# on a 2D domain for a 2D function $y=(y_1,y_2)$, a pressure $p$,
!# a dual velocity $\lambda$ and a dual pressure $\xi$. $u$ is the control
!# and $z$ a desired flow field.
!#
!# on a 2D domain for a 2D velocity function $y=(y_1,y_2)$, a pressure $p$,
!# a 2D dual velocity $\lambda=(\lambda_1,\lambda_2)$ and a dual pressure
!# $\xi$.
!#
!# A linear block system with 6 solution components is set up,
!# discretised and solved.
!#
!# </purpose>
!##############################################################################

program cc2doptc

  use main_program
  
  implicit none
  
  call main_optc()

end program
