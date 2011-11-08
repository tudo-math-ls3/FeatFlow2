!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic problem definitions for the cc2dmini_method2
!# solver. The basic structure and content of the different structures
!# are described here.
!# </purpose>
!##############################################################################

module cc2dminim2basic

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use boundarycondition
  use discretebc
  use linearsystemscalar
  use linearsystemblock
  
  use collection
    
  implicit none
  
  ! Maximum allowed level in this application; must be =9 for
  ! FEAT 1.x compatibility (still)!
  integer, parameter :: NNLEV = 9
  
!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! A system matrix for that specific level.
    type(t_matrixBlock) :: rmatrix

    ! Laplace matrix for that specific level.
    type(t_matrixScalar) :: rmatrixLaplace

    ! B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2

    ! A temporary vector for building the solution when assembling the
    ! matrix on lower levels.
    type(t_vectorBlock) :: rtempVector

    ! A variable describing the discrete boundary conditions fo the velocity
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for Nav.St. problem">

  type t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: NLMIN
    
    ! Maximum refinement level
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! A variable describing the analytic boundary conditions.
    type(t_boundaryConditions), pointer :: p_rboundaryConditions

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    type(t_problem_lvl), dimension(NNLEV) :: RlevelInfo
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection
    
    ! A param list that saves all parameters from the DAT/INI file(s).
    type(t_parlist) :: rparamList

  end type

!</typeblock>

!</types>

!******************************************************************************
! Documentation
!******************************************************************************

!******************************************************************************
! The problem structure t_problem
!******************************************************************************
!
! The problem structure t_problem collects all information of importance
! for the main CC2D solver module. It contains:
!  - Analytic definition of the boundary
!  - Analytic boundary conditions
!  - Information from the INI/DAT files
!  - Main solution and RHS vector
!  - Information for preconditioning in nonlinear loop
!  - For every level in the discretisation:
!    - Triangulation
!    - Block discretisation
!    - Matrices and
!    - Temporary vectors
! This problem structure is available only in the top-level routines of
! the CC2D module; it's not available in any callback routines.
! Parameters for callback routines are passed through the collection
! structure rcollection, which is part of the problem structure.
!
!
!******************************************************************************
! The collection structure t_problem%rcollection
!******************************************************************************
!
! The collection structure collects all information that is necessary to be
! passed to callback routines. This e.g. allows to pass matrices/vectors/
! constants from the main problem.
! This structure contains the following informaion, which is added by the
! initialisation routines to the collection:
!
! Global, level independent data:
!
! Name                  | Description
! ----------------------+------------------------------------------------------
! INI                   | t_paramlist object
!                       | Contains parameters from the DAT/INI files
!                       |
! NU                    | Reciprocal 1/RE of parameter RE from the DAT file
! NLMIN                 | Minimum level of the discretisation
! NLMAX                 | Maximum level of the discretisation
! ISTOKES               | =0: we discretise Navier Stokes,
!                       | =1: we discretise Stokes
! IUPWIND               | Type of stabilisation. 0=streamline diff, 1=upwind
! UPSAM                 | Stabilisation parameter
!
! On every level between NLMIN and NLMAX:
!
! Name                  | Description
! ----------------------+------------------------------------------------------
! LAPLACE               | Laplace matric
! SYSTEMMAT             | Nonlinear system matrix
! RTEMPVEC              | Temporary vector, compatible to matrix
!
!
! Global, level independent data, available during the nonlinear iteration:
!
! Name                  | Description
! ----------------------+------------------------------------------------------
! RHS                   | t_vectorBlock object
!                       | Current RHS vector on maximum level
!                       |
! SOLUTION              | t_vectorBlock object
!                       | Current solution vector on maximum level
!                       |
! ILVPROJECTION         | t_interlevelProjectionBlock structure.
!                       | Configures prolongation/restriction for multigrid
!                       | solver component.
!                       |
! RTEMPSCALAR           | t_vectorScalar object
!                       | Temporary vector
!                       |
! RTEMP2SCALAR          | t_vectorScalar object
!                       | Temporary vector
!                       |
! LINSOLVER             | t_linsol object
!                       | Configures the linear solver for preconditioning
!                       |
! OMEGAMIN              | Minimum damping parameter for correction in nonlinear
!                       | solver
! OMEGAMAX              | Maximum damping parameter for correction in nonlinear
!                       | solver
! IADAPTIVEMATRIX       | Whether to activate adaptive matrix construction
!                       | for coarse grid matrices
! DADMATTHRESHOLD       | Threshold parameter dor adaptive matrix construction
  
end module
