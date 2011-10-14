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

  private
  public :: t_perfconfig
  public :: pcfg_initPerfConfig

!<types>

!<typeblock>

  ! Global performance configuration
  type t_perfconfig

    ! Number of equations to be handled simultaneously
    integer :: NEQSIM    = 32

    ! Number of matrix entries to be handled simultaneously
    integer :: NASIM     = 32

    ! Number of edges to be handled simultaneously
    integer :: NEDGESIM  = 32

    ! Number of elements to be handled simultaneously
    integer :: NELEMSIM  = 128

    ! Number of patches to be handled simultaneously
    integer :: NPATCHSIM = 100

    ! Number of items to be handles simultaneously
    integer :: NITEMSIM  = 256

    ! OpenMP-Extension: the following settings are lower bounds which
    ! must be satisfied before OpenMP-parallelisation is activated

    ! Minimal number of equations
    !$ integer :: NEQMIN_OMP    = 1000

    ! Minimal number of matrix entries
    !$ integer :: NAMIN_OMP     = 1000

    ! Minimal number of edges
    !$ integer :: NEDGEMIN_OMP  = 1000

    ! Minimal number of elements
    !$ integer :: NELEMMIN_OMP  = 1000

    ! Minimal number of patches
    !$ integer :: NPATCHMIN_OMP = 1000

    ! Minimal number of items
    !$ integer :: NITEMMIN_OMP  = 1000

  end type

!</typeblock>

!</types>

contains

  !****************************************************************************

!<subroutine>

  subroutine pcfg_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<outpu>
  ! Performance configuration to be initialised
    type(t_perfconfig), intent(out) :: rperfconfig
!</input>
!</subroutine>
    
    rperfconfig%NEQSIM   = 32
    rperfconfig%NASIM    = 32
    rperfconfig%NEDGESIM = 32
    rperfconfig%NELEMSIM = 128
    rperfconfig%NPATCHSIM = 100
    rperfconfig%NITEMSIM  = 256

    ! OpenMP-Extension
    !$ rperfconfig%NEQMIN_OMP    = 1000
    !$ rperfconfig%NAMIN_OMP     = 1000
    !$ rperfconfig%NEDGEMIN_OMP  = 1000
    !$ rperfconfig%NELEMMIN_OMP  = 1000
    !$ rperfconfig%NPATCHMIN_OMP = 1000
    !$ rperfconfig%NITEMMIN_OMP  = 1000
    
  end subroutine pcfg_initPerfConfig

end module perfconfig
