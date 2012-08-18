!##############################################################################
!# ****************************************************************************
!# <name> CahnHilliard_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic problem structures for the heat conduction
!# problem.
!# </purpose>
!##############################################################################

module CahnHilliard_basic

  use fsystem
  use linearsolver
  use boundary
  use matrixfilters
  use vectorfilters
  use triangulation
  use spatialdiscretisation
  use sortstrategy
  use discretebc
  use discretefbc
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  use filtersupport
  
! Take care of this ...how to use adaptive time stepping?
  use timestepping
  use adaptivetimestep
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

  use collection
  use paramlist
    
  IMPLICIT NONE
  
!<types>

!<typeblock description="type block defining all information about one level">

  type t_CHproblem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation => NULL()

!MCai, add a rdiscretisation
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A scalar discretisation structure that specifies how to generate
    ! the mass matrix in the velocity FEM space. Take care.
    type(t_spatialDiscretisation) ::rdiscretisationLaplace
    type(t_spatialDiscretisation) ::rdiscretisationMass

    ! A template FEM matrix that defines the structure of A, B, C, D
    ! matrices. The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateFEM

   ! A-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixA

   ! B-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB

    ! C-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixC
    
    ! D-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixD

    ! System matrix. May change during the time iteration
    type(t_matrixBlock) :: rmatrix
    
    ! The mass matrix for phase variable (also chem potential)
    type(t_matrixScalar)  :: rmatrixMass

    ! The Laplace matrix for phase variable (also chem potential)
    type(t_matrixScalar) :: rmatrixLaplace

    ! The convection matrix (velocity \cdot grad \phi)
	type(t_matrixScalar) :: rmatrixConv

    ! We may need temp vector in nonlinear iteration.
    type(t_vectorBlock) :: rtempVector
    
    ! A variable describing the discrete boundary conditions.
    ! The variable points to NULL until the first boundary conditions
    ! are assembled.
    type(t_discreteBC), POINTER :: p_rdiscreteBC => NULL()
  
  end type
  
!</typeblock>

!<typeblock description="Configuration block for the time stepping.">

  type t_CHproblem_nonst

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Configuration block of the time stepping scheme.
    type(t_explicitTimeStepping)        :: rtimestepping
    
    ! Configuration block for the adaptive time stepping.
    type(t_adaptimeTimeStepping) :: radaptiveTimeStepping
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Number of current time step
    integer                             :: iiteration
    
    ! Maximum number of time steps
    integer                             :: niterations
    
    ! Start time
    real(DP)                            :: dtimemin
    
    ! Current time
    real(DP)                            :: dtime
    
    ! Maximum time
    real(DP)                            :: dtimemax
  
  end type

!</typeblock>

! Do we need adaptive time step? if yes, please add a t_CHproblem_explTimeStepping,
! similar to cc2d.

!<typeblock description="Application-specific type block for CahnHilliard problem">

  type t_CHproblem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: NLMIN
    
    ! Maximum refinement level
    integer :: NLMAX

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A RHS vector on the finest level used for solving linear systems
    type(t_vectorBlock) :: rrhs
    
    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), POINTER :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! A filter chain to filter vectors during the solution process
    type(t_filterChain), DIMENSION(1) :: RfilterChain

    integer :: itimedependence=1

    ! A parameter block for everything that controls the time dependence.
    type(t_CHproblem_nonst) :: rtimedependence

    ! An array of t_CHproblem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    type(t_CHproblem_lvl), DIMENSION(:), POINTER :: RlevelInfo

    ! A param list that saves all parameters from the DAT/INI file(s).
    type(t_parlist)                       :: rparamList

    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection
    
  end type

!</typeblock>

!</types>

end module