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
  
  ! Return the evaluation tag of the function
  integer, parameter, public :: UFCI_TASK_GET_EVAL_TAG          = 2

!</constantblock>

!<constantblock description="Callback function evaluation tags">

  ! Function needs element index list
  integer(I32), parameter, public :: UFCI_EVALTAG_ELEM_IDX      = 2_I32**0

  ! Function needs reference coordinates
  integer(I32), parameter, public :: UFCI_EVALTAG_REF_COORDS    = 2_I32**1
  
  ! Function needs real coordinates
  integer(I32), parameter, public :: UFCI_EVALTAG_REAL_COORDS   = 2_I32**2
  
  ! Function needs jacobian matrices
  integer(I32), parameter, public :: UFCI_EVALTAG_JAC           = 2_I32**3
  
  ! Function needs jacobian determinants
  integer(I32), parameter, public :: UFCI_EVALTAG_DETJ          = 2_I32**4

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
    
    ! ---------- ENTRIES FOR EVAL-TAG TASK ----------

    ! OUT: The evaluation tag of the functions. If nfunctions is > 1, then
    !      cevalTag is the OR-ed combination of all called functions.
    integer(I32) :: cevalTag = UFCI_EVALTAG_REAL_COORDS
    
    ! ---------- ENTRIES FOR EVALUATION TASK ----------

    ! IN: Total number of elements
    integer :: nelements = 0

    ! IN: Total number of points per element
    integer :: npointsPerElement = 0

    ! IN: A pointer to the triangulation structure on which the current
    !     evaluation is performed; is always given.
    type(t_triangulation), pointer :: p_rtria => null()
    
    ! IN: Element index list; may be null if the UFCI_EVALTAG_ELEM_IDX tag was
    !     not set.
    integer, dimension(:), pointer :: p_IelementIdx => null()

    ! IN: Points where to evaluate, given in reference coordinates; may be
    !     null if the UFCI_EVALTAG_REF_COORDS tag was not specified.
    real(DP), dimension(:,:,:), pointer :: p_DpointsRef => null()

    ! IN: Points where to evaluate, given in real coordinates; may be null if
    !     the UFCI_EVALTAG_REAL_COORDS tag was not specified.
    real(DP), dimension(:,:,:), pointer :: p_DpointsReal => null()

    ! IN: Jacobian matrices in the evaluation points; may be null if the
    !     UFCI_EVALTAG_JAC tag was not specified.
    real(DP), dimension(:,:,:), pointer :: p_Djac => null()

    ! IN: Jacobian determinants in the evaluation points; may be null if the
    !     UFCI_EVALTAG_DETJ tag was not specified.
    real(DP), dimension(:,:), pointer :: p_Ddetj => null()

    ! OUT: The values of the evaluated functions, where
    !      p_Dvalues(i,j,k) = function value of function k on element j in
    !                         cubature point i
    real(DP), dimension(:,:,:), pointer :: p_Dvalues => null()
   
  end type

!</typeblock>

!</types>

end module