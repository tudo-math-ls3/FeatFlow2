!##############################################################################
!# ****************************************************************************
!# <name> ufci </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the necessary structures and definitions which are
!# necessary to call the new 'unified function callback interface' defined in
!# 'intf_ufci.inc'.
!# </purpose>
!##############################################################################

module ufci

use fsystem
use triangulation

implicit none

private

public :: t_ufciData

!<constants>

!<constantblock description="General parameters">

  ! Maximum number of functions which may be evaluated simultaneously
  integer, parameter, public :: UFCI_MAX_SIM                    = 32

!</constantblock>

!<constantblock description="Task identifiers for callback function evaluation">

  ! Do nothing - this is just a dummy
  integer, parameter, public :: UFCI_TASK_NONE                  = 0

  ! Evaluate the function in a given set of points
  integer, parameter, public :: UFCI_TASK_EVALUATE              = 1

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Data structure for callback functions
  type t_ufciData

    ! IN: The task that is to be performed by the callback function; one of 
    !     the UFCI_TASK_XXXX constants which specifies what information is
    !     desired from the callback function by the caller.
    integer :: ctask = UFCI_TASK_NONE
    
    ! IN: Total number of functions which are to be evaluated.
    integer :: nfunctions = 0

    ! IN: user-defined ID of the functions which are to be called
    integer, dimension(UFCI_MAX_SIM) :: IfunctionIDs = 0
    
    ! ---------- ENTRIES FOR EVALUATION TASK ----------

    ! IN: The total number of points.
    integer :: npoints = 0

    ! IN: The points in which the functions are to be evaluated, given in real
    !     coordinates.
    real(DP), dimension(:,:), pointer :: p_Dpoints => null()

    ! OUT: The values of the evaluated functions, where
    !      p_Dvalues(i,j) = function value of function j in point i
    real(DP), dimension(:,:), pointer :: p_Dvalues => null()
   
  end type

!</typeblock>

!</types>

end module