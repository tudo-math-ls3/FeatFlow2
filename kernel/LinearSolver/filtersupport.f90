!##############################################################################
!# ****************************************************************************
!# <name> filtersupport </name>
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
!# The routine filter_applyFilterChainVec can be used to apply such a filter
!# chain to a vector. The routine stops when either the end of the array
!# is reached or if a filter is found that contains the tag FILTER_NOFILTER.
!#
!# A filter is always applied to the whole solution/defect vector. One filter
!# is for example the "apply-discrete-boundary-conditions" filter
!# that loops about all subvectors of a given vector to implement the
!# discrete boundary conditions. 
!#
!# The following routines can be found here:
!#
!# 1.) filter_applyFilterChainVec
!#     -> Applies a given filter chain onto a (block) vector.
!#
!# 2.) filter_applyFilterChainMat
!#     -> Applies a given filter chain onto a (block) matrix.
!# </purpose>
!##############################################################################

MODULE filtersupport

  USE fsystem
  USE linearsystemblock
  USE spatialdiscretisation
  USE matrixfilters
  USE vectorfilters
  
  IMPLICIT NONE

!<constants>

!<constantblock description="General constants concerning filters">

  ! A standard length for arrays holding a filter chain.
  INTEGER, PARAMETER :: FILTER_MAXFILTERS        =  32

!</constantblock>

!<constantblock description="Filter type flags for t\_filterChain%ifilterType">

  ! Undefined filter
  INTEGER, PARAMETER :: FILTER_NOFILTER          =  0

  ! Do-nothing filter; does nothing to a vector
  INTEGER, PARAMETER :: FILTER_DONOTHING         =  0

  ! Vector filter for imposing discrete boundary conditions of the 
  ! real boundary into a solution vector.
  INTEGER, PARAMETER :: FILTER_DISCBCSOLREAL     =  1

  ! Vector filter for imposing discrete boundary conditions of the 
  ! real boundary into a right-hand-side vector.
  INTEGER, PARAMETER :: FILTER_DISCBCRHSREAL     =  2

  ! Vector filter for imposing discrete boundary conditions of the 
  ! real boundary into a defect vector.
  INTEGER, PARAMETER :: FILTER_DISCBCDEFREAL     =  3

  ! Matrix filter for imposing discrete boundary conditions of the 
  ! real boundary into a matrix.
  INTEGER, PARAMETER :: FILTER_DISCBCMATREAL     =  4

  ! Vector filter for imposing discrete boundary conditions of 
  ! the fictitious boundary into a solution vector.
  INTEGER, PARAMETER :: FILTER_DISCBCSOLFICT     =  5

  ! Vector filter for imposing discrete boundary conditions of 
  ! the fictitious boundary into a right-hand-side vector.
  INTEGER, PARAMETER :: FILTER_DISCBCRHSFICT     =  6

  ! Vector filter for imposing discrete boundary conditions of 
  ! the fictitious boundary into a defect vector.
  INTEGER, PARAMETER :: FILTER_DISCBCDEFFICT     =  7

  ! Matrix filter for imposing discrete boundary conditions of the 
  ! fictitious boundary into a matrix.
  INTEGER, PARAMETER :: FILTER_DISCBCMATFICT     =  8

  ! Vector filter for bringing a subvector of a vector to the space $L^2_0$.
  INTEGER, PARAMETER :: FILTER_TOL20             =  9

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
    
    ! Information tag for the TOL20 filterif ifilterType=FILTER_TOL20:
    ! Number of the subvector that should be filtered to be in the
    ! space $L^2_0$.
    INTEGER                            :: itoL20component = 0
    
  END TYPE
  
  !</typeblock>

!</types>
  
CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE filter_applyFilterChainVec (rx, RfilterChain)

!<description>
  ! This routine applies a filter chain on a (block) vector rx. All filters in
  ! the chain are applied one after the other until the end of the array is 
  ! reached or the first filter is found with the filter tag FILTER_NOFILTER.
!</description>

!<input>
  
  ! The filter chain
  TYPE(t_filterChain), DIMENSION(:), INTENT(IN) :: RfilterChain
  
!</input>

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
  
    ! Choose the filter and apply it
    SELECT CASE (ifilterType)
    CASE (FILTER_NOFILTER)
      ! Cancel if we reached the last filter before reaching the end of the
      ! array
      EXIT
      
    CASE (FILTER_DISCBCSOLREAL)
      ! Impose Dirichlet boundary contitions into the solution vector rx
      CALL vecfil_discreteBCsol (rx)
    
    CASE (FILTER_DISCBCRHSREAL)
      ! Impose Dirichlet boundary contitions into the RHS vector rx
      CALL vecfil_discreteBCrhs (rx)
    
    CASE (FILTER_DISCBCDEFREAL)
      ! Impose Dirichlet boundary contitions into the defect vector rx
      CALL vecfil_discreteBCdef (rx)

    CASE (FILTER_DISCBCSOLFICT)
      ! Impose Dirichlet fictitious boundary contitions into the solution vector rx
      CALL vecfil_discreteFBCsol (rx)
    
    CASE (FILTER_DISCBCRHSFICT)
      ! Impose Dirichlet fictitious boundary contitions into the RHS vector rx
      CALL vecfil_discreteFBCrhs (rx)
    
    CASE (FILTER_DISCBCDEFFICT)
      ! Impose Dirichlet fictitious boundary contitions into the defect vector rx
      CALL vecfil_discreteFBCdef (rx)

    CASE (FILTER_TOL20)
      ! Bring the subvector itoL20component of rx to the space $L^2_0$:
      CALL vecfil_subvectorToL20 (rx,RfilterChain(i)%itoL20component)

    CASE DEFAULT
      PRINT *,'filter_applyFilterChainVec: Unknown filter.'
      EXIT
    END SELECT
  
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE filter_applyFilterChainMat (rmatrix, RfilterChain)

!<description>
  ! This routine applies a filter chain on a (block) matrix rmatrix. All 
  ! filters in the chain are applied one after the other until the end of 
  ! the array is reached or the first filter is found with the filter 
  ! tag FILTER_NOFILTER.
!</description>

!<input>
  
  ! The filter chain
  TYPE(t_filterChain), DIMENSION(:), INTENT(IN) :: RfilterChain
  
!</input>

!<inputoutput>
  
  ! The matrix where the filter will be applied to.
  ! This is also the result matrix.
  TYPE(t_matrixBlock), INTENT(INOUT)            :: rmatrix
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i, ifilterType
  
  ! Apply the filters to matrixx - one after the other
  
  DO i=LBOUND(RfilterChain,1),UBOUND(RfilterChain,1)
  
    ifilterType = RfilterChain(i)%ifilterType
  
    ! Choose the filter and apply it
    SELECT CASE (ifilterType)
    CASE (FILTER_NOFILTER)
      ! Cancel if we reached the last filter before reaching the end of the
      ! array
      EXIT
      
    CASE (FILTER_DISCBCSOLREAL,FILTER_DISCBCRHSREAL, &
          FILTER_DISCBCDEFREAL,FILTER_DISCBCMATREAL)
      ! Impose Dirichlet boundary contitions into the matrix rmatrix.
      ! The filter is the same for both, solution and defect filter,
      ! as the matrix modification is teh same (for now).
      CALL matfil_discreteBC (rmatrix)
    
    CASE (FILTER_DISCBCSOLFICT,FILTER_DISCBCRHSFICT, &
          FILTER_DISCBCDEFFICT,FILTER_DISCBCMATFICT)
      ! Impose Dirichlet fictitious boundary contitions into the matrix rmatrix.
      ! The filter is the same for both, solution and defect filter,
      ! as the matrix modification is teh same (for now).
      CALL matfil_discreteFBC (rmatrix)
    
    CASE DEFAULT
      PRINT *,'filter_applyFilterChainMat: Unknown filter.'
      EXIT
    
    END SELECT
  
  END DO
  
  END SUBROUTINE
  
END MODULE
