!##############################################################################
!# ****************************************************************************
!# <name> mprimitives </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of low-level mathmatic functions without
!# a big background of finite elements or similar, like linear transformation,
!# parabolic profile, etc.
!#
!# The routines here are usually declared as PURE as they should not have
!# any side effects!
!#
!# The following routines can be found here:
!#
!# 1.) mprim_getParabolicProfile
!#     -> Calculates the value of a parabolic profile along a line.
!#
!# </purpose>
!##############################################################################

MODULE mprimitives

  USE fsystem
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<function>
  
  PURE REAL(DP) FUNCTION mprim_getParabolicProfile (dpos,dlength,dmaxvalue)
  
!<description>
  ! Calculates the value of a parabolic profile along a line.
  ! dpos is a parameter value along a line of length dlength. The return value
  ! is the value of the profile and has its maximum dmax at position 0.5.
!</description>
  
!<input>
  
  ! Position on the line segment where to calculate the velocity
  REAL(DP), INTENT(IN) :: dpos
  
  ! Length of the line segment
  REAL(DP), INTENT(IN) :: dlength
  
  ! Maximum value of the profile 
  REAL(DP), INTENT(IN) :: dmaxvalue
  
!</input>

!<result>
  ! The value of the parabolic profile at position $dpos \in [0,dmax]$.
!</result>
  
    ! Result: Value of the parbolic profile on position dpos on the line segment
    mprim_getParabolicProfile = 4*dmaxvalue*dpos*(dlength-dpos)/(dlength*dlength)
    
  END FUNCTION

END MODULE
