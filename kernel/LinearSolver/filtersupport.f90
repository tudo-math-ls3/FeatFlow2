!##############################################################################
!# ****************************************************************************
!# <name> FilterSupport </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains filter routines that can be used in the solution
!# process of the problem (especially when solving linear subproblems).
!# For using filters in a solver, a ''filter chain'' must be build up.
!# This is a simple array of filter structures that describe the filter
!# that should be applied to a vector one after the other.
!#
!# To build a filter chain, the application must declare an array
!# of type t_filterChain (e.g. of length FILTER_MAXFILTERS) and set the
!# filer type identifier.
!# The routine filter_applyFilterChain can be used to apply such a filter
!# chain to a vector. The routine stops when either the end of the array
!# is reached or if a filter is found that contains the tag FILTER_NOFILTER.
!#
!# A filter is always applied to the whole solution vector. One filter is
!# for example the "apply-discrete-boundary-conditions" filter
!# that loops about all subvectors of a given vector to implement the
!# discrete boundary conditions. 
!# </purpose>
!##############################################################################

MODULE filtersupport

  USE fsystem
  USE linearsystemblock
  USE spatialdiscretisation
  USE vectorfilters
  
  IMPLICIT NONE

!<constants>

  !<constantblock description="General constants concerning filters">

  ! A standard length for arrays holding a filter chain.
  INTEGER, PARAMETER :: FILTER_MAXFILTERS        =  32

  !</constantblock>

  !<constantblock description="Filter type flags for t_filterChain%ifilterType">

  ! Undefined filter
  INTEGER, PARAMETER :: FILTER_NOFILTER          =  0

  ! Do-nothing filter; does nothing to a vector
  INTEGER, PARAMETER :: FILTER_DONOTHING         =  0

  ! Vector filter for imposing discrete boundary conditions of the 
  ! real boundary into a solution vector.
  INTEGER, PARAMETER :: FILTER_DISCBCSOLREAL     =  1

  ! Vector filter for imposing discrete boundary conditions of the 
  ! real boundary into a defect vector.
  INTEGER, PARAMETER :: FILTER_DISCBCDEFREAL     =  2

  ! Vector filter for imposing discrete boundary conditions of 
  ! the fictitious boundary into a solution vector.
  INTEGER, PARAMETER :: FILTER_DISCBCSOLFICT     =  3

  ! Vector filter for imposing discrete boundary conditions of 
  ! the fictitious boundary into a defect vector.
  INTEGER, PARAMETER :: FILTER_DISCBCDEFFICT     =  4

  !</constantblock>
  
!</constants>

!<types>
  
  !<typeblock>
  
  ! A structure for using filters in a solver. This structure contains all
  ! information that are necessary to apply a filter (level independent!).
  
  TYPE t_filterChain
    
    ! Tag that identifies the type of filter to be applied to a vector.
    ! One of the FILTER_XXXX constants
    INTEGER                            :: ifilterType = FILTER_NOFILTER
    
  END TYPE
  
  !</typeblock>

!</types>
  
CONTAINS

!<subroutine>

  SUBROUTINE filter_applyFilterChain (rx, RfilterChain)

  !<description>
  
  ! This routine applies a filter chain on a vector rx. All filters in the
  ! chain are applied one after the other until the end of the array is reached
  ! or the first filter is found with the filter tag FILTER_NOFILTER.
  
  !</description>

  !<input>
  
  ! The filter chain
  TYPE(t_filterChain), DIMENSION(:), INTENT(IN) :: RfilterChain

  !<inputoutput>
  
  ! The vector where the filter will be applied to.
  ! This is also the result vector.
  TYPE(t_vectorBlock), INTENT(INOUT)            :: rx
  
  !</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i, ifilterType
  
  ! Apply the filters to rx - one after the other
  
  DO i=LBOUND(RfilterChain,1),UBOUND(RfilterChain,1)
  
    ifilterType = RfilterChain(i)%ifilterType
  
    ! Cancel if we reached the last filter before reaching the end of the
    ! array
    IF (ifilterType .EQ. FILTER_NOFILTER) EXIT
    
    ! Choose the filter and apply it
    SELECT CASE (ifilterType)
    CASE (FILTER_NOFILTER)
      ! Cancel if we reached the last filter before reaching the end of the
      ! array
      EXIT
    CASE (FILTER_DISCBCDEFREAL)
      ! Impose Dirichlet boundary contitions into the defect vector rx
      CALL vecfil_imposeDiscreteDefectBC (rx)
    END SELECT
  
    ! Apply the filter
    PRINT *,'to be implemented'
  
  END DO
  
  END SUBROUTINE
  
END MODULE
