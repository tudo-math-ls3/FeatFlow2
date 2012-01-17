!##############################################################################
!# ****************************************************************************
!# <name> fictitiousboundary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module defines the basic structure to handle fictitious boundary
!# objects.
!#
!# A fictitious boundary object is a general, analytically described part of
!# the underlying domain. To identify it, the structure t_fictBoundaryRegion
!# is defined, which contains basically some tags like a string for
!# the name of the fictitious boundary region. The analytic description
!# itself must be given by the application.
!# </purpose>
!##############################################################################

module fictitiousboundary

!$use omp_lib
  use fsystem

  implicit none
  
  private

!<types>

!<typeblock>

  ! Structure for defining a fictitious boundary region in the domain.
  ! The content is more or less application specific. The application must
  ! fill it with as much data as necessary to be able to identify a fictitious
  ! boundary region when discretising the associated boundary conditions.
  
  type t_fictBoundaryRegion
  
    ! Name of the region. Application specific.
    character(LEN=SYS_NAMELEN) :: sname = ''
    
    ! Integer identifier for the region. Application specific
    integer :: iidentifier = 0
  
    ! Integer value containing some flags that specify additional
    ! properties of the region. Application specific.
    integer(I32) :: iflags = 0
    
    ! Integer tag. Application specific
    integer :: itag = 0
    
    ! Double precision tag. Application specific
    real(DP) :: dtag = 0.0_DP
    
  end type
  
  public :: t_fictBoundaryRegion
  
!</typeblock>

! </types>

end module

