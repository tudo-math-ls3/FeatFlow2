!##############################################################################
!# ****************************************************************************
!# <name> perfconfig </name>
!# ****************************************************************************
!#
!# <purpose>
!#   This module provides basic structures for performance configurations
!# </purpose>
!##############################################################################

module perfconfig

  implicit none

  public :: t_perfconfig
  private

!<types>

!<typeblock>

  ! Global performance configuration
  type t_perfconfig

    ! Number of equations to be handled simultaneously
    integer :: NEQSIM   = 32

    ! Number of matrix entries to be handled simultaneously
    integer :: NASIM    = 32

    ! Number of edges to be handled simultaneously
    integer :: NEDGESIM = 32

    ! Number of elements to be handled simultaneously
    integer :: NELEMSIM = 128

    ! Number of patches to be handled simultaneously
    integer :: NPATCHSIM = 100

    ! Number of items to be handles simultaneously
    integer :: NITEMSIM  = 256

    ! OpenMP-Extension: the following settings are lower bounds which
    ! must be satisfied before OpenMP-parallelisation is activated

    ! Minimal number of equations
    !$ integer :: NEQMIN_OMP   = 1000

    ! Minimal number of matrix entries
    !$ integer :: NAMIN_OMP    = 1000

    ! Minimal number of edges
    !$ integer :: NEDGEMIN_OMP = 1000

    ! Minimal number of elements
    !$ integer :: NELEMMIN_OMP = 1000

    ! Minimal number of patches
    !$ integer :: NPATCHMIN_OMP = 1000

    ! Minimal number of items
    !$ integer :: NITEMMIN_OMP  = 1000

  end type

!</typeblock>

!</types>

end module perfconfig
