!##############################################################################
!# ****************************************************************************
!# <name> x </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains
!# </purpose>
!##############################################################################

module x

!$use omp_lib
  use fsystem
  
  implicit none

  private
  public :: x

!<types>
  
  !<typeblock>
  
  type x
    !...
  end type
  
  !</typeblock>

!</types>
  
  ! ***************************************************************************
  
!<constants>

  !<constantblock description="x">

  ! logical value 'true'
  integer, parameter, public :: YES = 0

  !</constantblock>
  
!</constants>
 
  ! ***************************************************************************


end module
