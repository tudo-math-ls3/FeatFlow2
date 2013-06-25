!##############################################################################
!# ****************************************************************************
!# <name> LS_Poisson_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the poisson problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# 
!# The function corresponds to a specific interface defined in 'intf_xxxx.inc'
!# files.
!#
!# 1.) getBoundaryValues_2D
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# </purpose>
!##############################################################################

module LS_Poisson_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use spatialdiscretisation
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use genoutput
  use triangulation
  use derivatives
  use linearsystemscalar
  use linearsystemblock
  use element
  use paramlist
  
  implicit none

contains

  ! ***************************************************************************


!<subroutine>

  subroutine getBoundaryValues_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...)
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(in)                                         :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                                         :: cinfoNeeded
  
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  integer, intent(in)                                          :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  real(DP), intent(in)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>

 ! Local variable
    real(DP) :: n
    ! parameter
    type(t_parlist) :: rparams 
  
!</subroutine>

!   ! To get the X/Y-coordinates of the boundary point, use:
!   REAL(DP) :: dx,dy
!   CALL boundary_getCoords(rdiscretisation%p_rboundary, &
!         rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
          
  
!   !     This is the 1st Edge
!     if ((dwhere.ge.0.0_DP).and.(dwhere.le.1.0_DP)) then
!     Dvalues(1) = 8.0_DP*dx*dx
!     end if
! !     2nd Edge
!     if ((dwhere.ge.1.0_DP).and.(dwhere.le.2.0_DP)) then
!     Dvalues(1) = 8.0_DP*(1.0_DP+dy*dy)
!     end if
! !     3rd Edge
!     if ((dwhere.ge.2.0_DP).and.(dwhere.le.3.0_DP)) then
!     Dvalues(1) = 8.0_DP*(1.0_DP+dx*dx)
!     end if
! !     4th Edge
!     if ((dwhere.ge.3.0_DP).and.(dwhere.le.4.0_DP)) then
!     Dvalues(1) = 8.0_DP*(dy*dy)
!     end if

  Dvalues(1) = 0.0_DP

  end subroutine
 

end module
