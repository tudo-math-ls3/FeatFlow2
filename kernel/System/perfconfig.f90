!##############################################################################
!# ****************************************************************************
!# <name> perfconfig </name>
!# ****************************************************************************
!#
!# <purpose>
!#   This module provides basic structures for performance configurations
!#
!# The following routines can be found in this module:
!#
!# 1.) pcfg_initFromParameterList
!#     -> Initialises a performance configuration by the values of a
!#        given parameter list
!#
!# 2.) pcfg_readPerfConfig
!#     -> Reads a performance configuration from file
!#
!# 3.) pcfg_writePerfConfig
!#     -> Writes a performance configuration to file
!# </purpose>
!##############################################################################

module perfconfig

!$use omp_lib
  use fsystem
  use paramlist
  
  implicit none

  private
  public :: t_perfconfig
  public :: pcfg_initPerfConfig
  public :: pcfg_initFromParameterList
  public :: pcfg_readPerfConfig
  public :: pcfg_writePerfConfig

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

!<output>
    ! Performance configuration to be initialised
    type(t_perfconfig), intent(out) :: rperfconfig
!</output>
!</subroutine>
    
    rperfconfig%NEQSIM    = 32
    rperfconfig%NASIM     = 32
    rperfconfig%NEDGESIM  = 32
    rperfconfig%NELEMSIM  = 128
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

  !****************************************************************************

!<subroutine>

  subroutine pcfg_initFromParameterlist(rparlist, ssectionName, rperfconfig)

!<description>
    ! This subroutine initialises a performance configuration by the
    ! values from a given parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in)    :: rparlist

    ! Section name of the parameter list
    character(LEN=*), intent(in)   :: ssectionName
!</input>

!<output>
  ! Performance configuration to be initialised
    type(t_perfconfig), intent(out) :: rperfconfig
!</output>
!</subroutine>

    ! Initialise xxxSIM settings
    call parlst_getvalue_int(rparlist, ssectionName,&
        "NEQSIM", rperfconfig%NEQSIM, rperfconfig%NEQSIM)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "NASIM", rperfconfig%NASIM, rperfconfig%NASIM)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "NEDGESIM", rperfconfig%NEDGESIM, rperfconfig%NEDGESIM)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "NELEMSIM", rperfconfig%NELEMSIM, rperfconfig%NELEMSIM)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "NPATCHSIM", rperfconfig%NPATCHSIM, rperfconfig%NPATCHSIM)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "NITEMSIM", rperfconfig%NITEMSIM, rperfconfig%NITEMSIM)

    ! OpenMP-Extension: xxxMIN_OMP
    !$ call parlst_getvalue_int(rparlist, ssectionName,&
    !$    "NEQMIN_OMP", rperfconfig%NEQMIN_OMP, rperfconfig%NEQMIN_OMP)
    !$ call parlst_getvalue_int(rparlist, ssectionName,&
    !$    "NAMIN_OMP", rperfconfig%NAMIN_OMP, rperfconfig%NAMIN_OMP)
    !$ call parlst_getvalue_int(rparlist, ssectionName,&
    !$    "NEDGEMIN_OMP", rperfconfig%NEDGEMIN_OMP, rperfconfig%NEDGEMIN_OMP)
    !$ call parlst_getvalue_int(rparlist, ssectionName,&
    !$    "NELEMMIN_OMP", rperfconfig%NELEMMIN_OMP, rperfconfig%NELEMMIN_OMP)
    !$ call parlst_getvalue_int(rparlist, ssectionName,&
    !$    "NPATCHMIN_OMP", rperfconfig%NPATCHMIN_OMP, rperfconfig%NPATCHMIN_OMP)
    !$ call parlst_getvalue_int(rparlist, ssectionName,&
    !$    "NITEMMIN_OMP", rperfconfig%NITEMMIN_OMP, rperfconfig%NITEMMIN_OMP)
    
  end subroutine pcfg_initFromParameterlist

  !****************************************************************************

!<subroutine>

  subroutine pcfg_readPerfConfig(sfilename, ssectionName, rperfconfig)

!<description>
  ! This routine reads a performance configuration from an INI-file
!</description>

!<input>
    ! File name from which the configuration is read
    character(len=*), intent(in) :: sfilename

    ! Section name of the INI-file
    character(LEN=*), intent(in)   :: ssectionName
!</input>

!<output>
    ! Performance configuration to be initialised
    type(t_perfconfig), intent(out) :: rperfconfig
!</output>
!</subroutine>

    ! local variables
    type(t_parlist) :: rparlist

    ! Read parameterlist from file
    call parlst_init(rparlist)
    call parlst_readfromfile(rparlist, sfilename)
    
    ! Initialise performance configuration
    call pcfg_initFromParameterList(rparlist, ssectionName, rperfconfig)

    ! Release parameter list
    call parlst_done(rparlist)

  end subroutine pcfg_readPerfConfig

  !****************************************************************************

!<subroutine>

  subroutine pcfg_writePerfConfig(rperfconfig, sfilename, ssectionName, cflag)

!<description>
  ! This routine writes a performance configuration to an INI-file
!</description>

!<input>
    ! Performance configuration to be initialised
    type(t_perfconfig), intent(in) :: rperfconfig

    ! File name from which the configuration is read
    character(len=*), intent(in) :: sfilename

    ! Section name of the INI-file
    character(LEN=*), intent(in)   :: ssectionName

    ! mode: SYS_APPEND or SYS_REPLACE
    integer, intent(in) :: cflag
!</input>
!</subroutine>

    ! local variables
    type(t_parlist) :: rparlist

    ! Initialise parameter list and add section
    call parlst_init(rparlist)
    call parlst_addsection(rparlist, ssectionname)

    ! Add xxxSIM constants
    call parlst_addvalue(rparlist, ssectionname,&
                         "NEQSIM", sys_siL(rperfconfig%NEQSIM,10))
    call parlst_addvalue(rparlist, ssectionname,&
                         "NASIM", sys_siL(rperfconfig%NASIM,10))
    call parlst_addvalue(rparlist, ssectionname,&
                         "NEDGESIM", sys_siL(rperfconfig%NEDGESIM,10))
    call parlst_addvalue(rparlist, ssectionname,&
                         "NELEMSIM", sys_siL(rperfconfig%NELEMSIM,10))
    call parlst_addvalue(rparlist, ssectionname,&
                         "NPATCHSIM", sys_siL(rperfconfig%NPATCHSIM,10))
    call parlst_addvalue(rparlist, ssectionname,&
                         "NITEMSIM", sys_siL(rperfconfig%NITEMSIM,10))

    ! OpenMP-Extension: xxxMIN_OMP constants
    !$ call parlst_addvalue(rparlist, ssectionname,&
    !$                      "NEQMIN_OMP", sys_siL(rperfconfig%NEQMIN_OMP,10))
    !$ call parlst_addvalue(rparlist, ssectionname,&
    !$                      "NAMIN_OMP", sys_siL(rperfconfig%NAMIN_OMP,10))
    !$ call parlst_addvalue(rparlist, ssectionname,&
    !$                      "NEDGEMIN_OMP", sys_siL(rperfconfig%NEDGEMIN_OMP,10))
    !$ call parlst_addvalue(rparlist, ssectionname,&
    !$                      "NELEMMIN_OMP", sys_siL(rperfconfig%NELEMMIN_OMP,10))
    !$ call parlst_addvalue(rparlist, ssectionname,&
    !$                      "NPATCHMIN_OMP", sys_siL(rperfconfig%NPATCHMIN_OMP,10))
    !$ call parlst_addvalue(rparlist, ssectionname,&
    !$                      "NITEMMIN_OMP", sys_siL(rperfconfig%NITEMMIN_OMP,10))

    ! Dump parameter list to file
    call parlst_dumptofile(rparlist, sfilename, cflag)

    ! Release parameter list
    call parlst_done(rparlist)

  end subroutine pcfg_writePerfConfig

end module perfconfig
