!##############################################################################
!# ****************************************************************************
!# <name> vanka </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the implementations of the VANKA preconditioner.
!# These are more or less auxiliary routines called by the VANKA
!# preconditioner in the linearsolver.f90 solver library.
!#
!# The following routines can be found here:
!#
!#  1.) vanka_initConformal
!#      -> Initialise the VANKA for conformal discretisations.
!#
!#  2.) vanka_conformal
!#      -> Perform one step of the VANKA solver conformal discretisations. Calls a
!#         specialised sub-VANKA-variant to do the work.
!#
!#  3.) vanka_doneConformal
!#      -> Clean up the VANKA for conformal discretisations.
!#
!# 4.) vanka_initPerfConfig
!#      -> Initialises the global performance configuration
!#
!# The following list of routines are all used internally. There is no need to
!# call them directly.
!#
!#  4.) vanka_initGeneralVanka
!#      -> Initialise the general Vanka solver
!#
!#  5.) vanka_general
!#      -> Perform one step of the full VANKA solver for general block systems.
!#
!#  6.) vanka_doneGeneralVanka
!#      -> Clean up the general VANKA solver
!#
!#  7.) vanka_init2DSPQ1TQ0simple
!#      -> Initialise specialised VANKA solver for 2D saddle point problems
!#         with <tex>$\tilde Q_1/Q_0$</tex> discretisation.
!#         Deprecated, ist not used.
!#
!#  8.) vanka_init2DNavierStokes
!#      -> Initialise the VANKA solver for the problem class
!#         '2D Navier-Stokes equation'.
!#
!#  9.) vanka_2DNavierStokes
!#      -> Apply the VANKA solver for the problem class
!#         '2D Navier-Stokes equation'.
!#
!# 10.) vanka_2DSPQ1TQ0simple
!#      -> Perform one step of the specialised VANKA solver for 2D saddle point
!#         problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation
!#
!# 11.) vanka_2DSPQ1TQ0simpleConf
!#      -> Perform one step of the specialised VANKA solver for 2D saddle point
!#         problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation.
!#         Applies VANKA only to a subset of all elements in the domain.
!#
!# 12.) vanka_2DSPQ1TQ0simpleCoupConf
!#      -> Perform one step of the specialised VANKA solver for 2D saddle point
!#         problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation.
!#         Applies VANKA only to a subset of all elements in the domain.
!#         This variant can handle fully coupled matrices.
!#
!# 13.) vanka_2DSPQ1TQ0fullConf
!#      -> Perform one step of the specialised 'full' VANKA solver for 2D saddle
!#         point problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation.
!#         Applies VANKA only to a subset of all elements in the domain.
!#
!# 14.) vanka_2DSPQ1TQ0fullCoupConf
!#      -> Perform one step of the specialised 'full' VANKA solver for 2D saddle
!#         point problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation.
!#         Applies VANKA only to a subset of all elements in the domain.
!#         This variant can handle fully coupled matrices.
!#
!# 15.) vanka_2DSPQ2QP1simple
!#      -> Perform one step of the specialised VANKA solver for 2D saddle point
!#         problems with <tex>$Q_2/QP_1$</tex> discretisation. Diagonal VANKA approach.
!#
!# 16.) vanka_2DSPQ2QP1full
!#      -> Perform one step of the specialised VANKA solver for 2D saddle point
!#         problems with <tex>$Q_2/QP_1$</tex> discretisation. Full VANKA approach.
!#
!# 17.) vanka_2DSPQ2QP1simpleConf
!#      -> Perform one step of the specialised VANKA solver for 2D saddle point
!#         problems with <tex>$Q_2/QP_1$</tex> discretisation. Diagonal VANKA approach.
!#         Applies VANKA only to a subset of all elements in the domain.
!#
!# 18.) vanka_2DSPQ2QP1fullConf
!#      -> Perform one step of the specialised VANKA solver for 2D saddle point
!#         problems with <tex>$Q_2/QP_1$</tex> discretisation. Full VANKA approach.
!#         Applies VANKA only to a subset of all elements in the domain.
!#
!# 19.) vanka_init2DNavierStokesOptC
!#      -> Initialise the VANKA solver for 2D Navier-Stokes optimal control
!#         problems. Specialised $\tilde Q1/Q0$ version, full VANKA approach.
!#
!# 20.) vanka_2DNSSOCQ1TQ0fullCoupConf
!#      -> Perform one step of the VANKA solver for 2D Navier-Stokes optimal
!#         control problems. Specialised $\tilde Q1/Q0$ version, full VANKA approach.
!#
!# 21.) vanka_init3DNavierStokes
!#      -> Initialise the VANKA solver for the problem class
!#         '3D Navier-Stokes equation'.
!#
!# 22.) vanka_3DNavierStokes
!#      -> Apply the VANKA solver for the problem class
!#         '3D Navier-Stokes equation'.
!#
!# 23.) vanka_3DSPQ1TQ0simple
!#      -> Perform one step of the specialised VANKA solver for 3D saddle point
!#         problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation
!#
!# 24.) vanka_3DSPQ1TQ0simpleConf
!#      -> Perform one step of the specialised VANKA solver for 3D saddle point
!#         problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation.
!#         Applies VANKA only to a subset of all elements in the domain.
!#
!# 25.) vanka_3DSPQ1TQ0fullConf
!#      -> Perform one step of the specialised 'full' VANKA solver for 3D saddle
!#         point problems with <tex>$\tilde Q_1/Q_0$</tex> discretisation.
!#         Applies VANKA only to a subset of all elements in the domain.
!#
!# 26.) vanka_2DSPQ2QP1simpleCoupConf
!#      -> VANKA for 2D Navier-Stokes with completely different A11/A12/A21/A22,
!#         Q2/QP1 discretisation.
!#
!#  History
!# ---------
!# Originally, the following VANKA variants were implemented:
!#  - vanka_general
!#  - vanka_init2DSPQ1TQ0simple
!# vanka_general was designed to work with everything. vanka_init2DSPQ1TQ0simple was
!# designed to work with the 2D Navier Stokes problem, uniform <tex>$\tilde Q_1/Q_0$</tex>
!# discretisation. It used the 'diagonal' VANKA approach.
!#
!# Over the time, there was a 'full' VANKA variant added as well as variants that
!# were designed to work with <tex>$Q_2/QP_1$</tex> discretisations. However, each variant
!# had to be chosen carefully, as all variants were rather specialised to a
!# special uniform discretisation.
!#
!# To work with more general configurations, a VANKA wrapper was introduced.
!# This replaces the needs for carefully choosing the correct VANKA variant
!# and can be seen as 'more general VANKA', as it is also designed to work
!# with conformal discretisations.
!#
!# For this purpose, all VANKA variants were modified to work with element
!# lists, e.g.
!#   vanka_2DSPQ1TQ0simple -> vanka_2DSPQ1TQ0simpleConf
!#   vanka_2DSPQ2QP1simple -> vanka_2DSPQ2QP1simpleConf
!#   ...
!# For the Navier-Stokes problem, the corresponding wrapper function is
!#  - vanka_initConformal / vanka_conformal
!# These will check the discretisation and call the correct
!# vanka_xxxxConf subroutines automatically.
!#
!# To support other types of problems, the user should proceed as follows:
!#  - Write a VANKA variant in the '....Conf' style that accepts an element
!#    list, which specifies a subset of the domain where to apply VANKA
!#    (like vanka_2DSPQ1TQ0simpleConf).
!#  - Add a VANKA subtype constant and a VANKA problem class constant
!#    to the VANKATP_xxxx and VANKACP_xxxx list if necessary.
!#  - Introduce new VANKA structures as necessary and add references
!#    to that to the t_vanka structure.
!#  - Modify the vanka_initConformal / vanka_conformal routines so that
!#    they initialise and call the new VANKA variant correctly.
!# </purpose>
!##############################################################################

module vanka

!$use omp_lib
  use dofmapping
  use element
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use perfconfig
  use quicksolver
  use spatialdiscretisation
  use storage

  use vanka_navst2d
  use vanka_bouss2d
  use vanka_optcontrol

  implicit none

  private

  public :: t_vanka
  public :: t_vankaGeneral
  public :: t_vankaPointer2DNavSt
  public :: t_matrixPointer79Vanka
  public :: t_vankaPointer2DNavStOptC
  public :: t_vankaPointer3DNavSt
  public :: vanka_initPerfConfig
  public :: vanka_initConformal
  public :: vanka_doneConformal
  public :: vanka_conformal

  ! Public entities imported from vanka_navst2d.f90
  public :: VANKATP_NAVST2D_DIAG
  public :: VANKATP_NAVST2D_FULL

  ! Public entities imported from vanka_bouss2d.f90
  public :: VANKATP_BOUSS2D_DIAG
  public :: VANKATP_BOUSS2D_FULL

  ! Public entities imported from vanka_optcontrol.f90
  public :: VANKATP_NAVSTOPTC2D_DIAG
  public :: VANKATP_NAVSTOPTC2D_FULL

!<constants>

!<constantblock description="Identifiers for the different problem classes that can be handled by VANKA">

  ! Dummy, Vanka not initialised.
  integer, parameter, public :: VANKAPC_NONE           = -1

  ! General VANKA
  integer, parameter, public :: VANKAPC_GENERAL        = 0

  ! 2D Navier Stokes problem
  integer, parameter, public :: VANKAPC_2DNAVIERSTOKES = 1

  ! 2D Navier Stokes optimal control problem
  integer, parameter, public :: VANKAPC_2DNAVIERSTOKESOPTC = 2

  ! 3D Navier Stokes problem
  integer, parameter, public :: VANKAPC_3DNAVIERSTOKES = 3

  ! --------------------- NEW IMPLEMENTATION ---------------------

  ! 2D Navier-Stokes problem
  integer, parameter, public :: VANKAPC_NAVIERSTOKES2D = 10

  ! 2D Boussinesq problem
  integer, parameter, public :: VANKAPC_BOUSSINESQ2D = 11

  ! 2D Navier-Stokes optimal control problem
  integer, parameter, public :: VANKAPC_NAVIERSTOKESOPTC2D = 12

!</constantblock>


!<constantblock description="Identifiers for the different VANKA subtypes. Which one is supported depends on the problem class.">

  ! Standard VANKA, most suitable for the corresponding situation
  integer, parameter, public :: VANKATP_STANDARD  = 0

  ! Diagonal-type VANKA
  integer, parameter, public :: VANKATP_DIAGONAL  = 0

  ! 'Full' VANKA
  integer, parameter, public :: VANKATP_FULL      = 1

  ! 'Full' VANKA for optimal control problems, primal equation processing
  integer, parameter, public :: VANKATP_FULLOPTC_PRIMAL = 2

  ! 'Full' VANKA for optimal control problems, dual equation processing
  integer, parameter, public :: VANKATP_FULLOPTC_DUAL   = 3

  ! 'Diagonal' VANKA for optimal control problems
  integer, parameter, public :: VANKATP_DIAGOPTC        = 4

  ! Diagonal-type VANKA, solution based
  integer, parameter, public :: VANKATP_DIAGONAL_SOLBASED = 5

!</constantblock>

!</constants>


!<types>

!<typeblock>

  ! A structure that accepts a pointer to the column/row/data arrays
  ! of a structure-7/structure 9 matrix. This is usually passed to
  ! the VANKA preconditioner(s) to specify the matrices to handle.


  type t_matrixPointer79Vanka
    ! Is set to FALSE if the matrix does not exist/is empty.
    ! In this case, the pointers below are undefined!
    logical :: bexists

    ! TRUE if the matrix is saved transposed
    logical :: btransposed

    ! The scaling factor of the matrix; from the matrix structure.
    real(DP) :: dscaleFactor

    ! Pointer to the data - currently only double precision
    real(DP), dimension(:), pointer :: p_DA

    ! Pointer to the column structure
    integer, dimension(:), pointer :: p_Kcol

    ! Pointer to the row structure
    integer, dimension(:), pointer :: p_Kld
  end type

!</typeblock>

!<typeblock>

  ! A structure for saving precalculated information for the general VANKA.
  ! This is initialised by vanka_initGeneralVanka and released by
  ! vanka_doneGeneralVanka.

  type t_vankaGeneral

    ! Number of blocks in the global matrix
    integer :: nblocks

    ! Pointer to the block matrix
    type(t_matrixBlock), pointer :: p_rmatrix

    ! Pointers to t_matrixPointer79Vanka structures specifying
    ! the submatrices and their properties.
    type(t_matrixPointer79Vanka), dimension(:,:),pointer :: p_Rmatrices

    ! Maximum number of local DOF`s.
    integer :: nmaxLocalDOFs

    ! Total number of local DOF`s
    integer :: ndofsPerElement

    ! Number of local DOF`s in the element distributions of all blocks.
    ! Note that this VANKA supports only uniform discretisations, so
    ! each entry corresponds to one block in the solution vector.
    integer, dimension(:), pointer :: p_InDofsLocal => null()

    ! Offset indices of the blocks in the solution vector. IblockOffset(i)
    ! points to the beginning of the i-th block of the solution vector.
    integer, dimension(:), pointer :: p_IblockOffset => null()

    ! Temporary array that saves the DOF`s that are in processing when
    ! looping over an element set.
    ! DIMENSION(nmaxLocalDOFs,NELEMSIM,nblocks)
    integer, dimension(:,:,:), pointer :: p_IelementDOFs => null()

  end type

!</typeblock>

!<typeblock>

  ! A structure that saves matrix pointers for the 2D-Navier-Stokes
  ! VANKA method for Navier-Stokes systems.

  type t_vankaPointer2DNavSt
    ! Pointer to the column structure of the velocity matrix A11 and A22
    ! A11 and A22 must have the same structure.
    integer, dimension(:), pointer :: p_KcolA => null()

    ! Pointer to the row structure of the velocity matrix A11 and A22.
    ! They must have the same structure.
    integer, dimension(:), pointer :: p_KldA => null()

    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KdiagonalA => null()

    ! Pointer to the matrix entries of the velocity matrix A11
    real(DP), dimension(:), pointer :: p_DA => null()

    ! Pointer to the matrix entries of the velocity matrix A22 or NULL
    ! A11=A22.
    real(DP), dimension(:), pointer :: p_DA22 => null()

    ! Pointer to the column structure of the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KcolA12 => null()

    ! Pointer to the row structure of the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KldA12 => null()

    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KdiagonalA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A12 or NULL
    ! if not present
    real(DP), dimension(:), pointer :: p_DA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A21 or NULL
    ! if not present
    real(DP), dimension(:), pointer :: p_DA21 => null()

    ! Pointer to the column structure of the B/D-matrices.
    integer, dimension(:), pointer :: p_KcolB => null()

    ! Pointer to the row structure of the B/D-matrices
    integer, dimension(:), pointer :: p_KldB => null()

    ! Pointer to the entries of the B1-matrix
    real(DP), dimension(:), pointer :: p_DB1 => null()

    ! Pointer to the entries of the B2-matrix
    real(DP), dimension(:), pointer :: p_DB2 => null()

    ! Pointer to the entries of the D1-matrix
    real(DP), dimension(:), pointer :: p_DD1 => null()

    ! Pointer to the entries of the D2-matrix
    real(DP), dimension(:), pointer :: p_DD2 => null()

    ! Pointer to the matrix entries of the pressure identity matrix A33
    ! (if it exists).
    real(DP), dimension(:), pointer :: p_DA33 => null()

    ! Pointer to diagonal entries of A33
    integer, dimension(:), pointer :: p_KdiagonalA33 => null()

    ! Spatial discretisation structure for X-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrU => null()

    ! Spatial discretisation structure for Y-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrV => null()

    ! Spatial discretisation structure for pressure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrP => null()

    ! Multiplication factors for the submatrices; taken from the system matrix.
    ! (-> Not used in the current implementation! Although it is easy to include
    ! that into VANKA, some further speed analysis has to be done to make
    ! sure there is not too much speed impact when using these!)
    real(DP), dimension(3,3) :: Dmultipliers

    ! A temporary vector for solution based VANKA variants.
    ! Undefined if not needed.
    type(t_vectorBlock) :: rtempVector

  end type

!</typeblock>

!<typeblock>

  ! A structure that saves matrix pointers for the 2D-Navier-Stokes
  ! VANKA method for Navier-Stokes optimal control systems.

  type t_vankaPointer2DNavStOptC
    ! Pointer to the column structure of the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KcolA11 => null()

    ! Pointer to the row structure of the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KldA11 => null()

    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KdiagonalA11 => null()

    ! Pointer to the matrix entries of the velocity matrix A11
    real(DP), dimension(:), pointer :: p_DA11 => null()

    ! Pointer to the matrix entries of the velocity matrix A22 or NULL
    ! if not present
    real(DP), dimension(:), pointer :: p_DA22 => null()

    ! Pointer to the column structure of the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KcolA12 => null()

    ! Pointer to the row structure of the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KldA12 => null()

    ! Pointer to diagonal entries in the velocity matrix A11 and A22
    integer, dimension(:), pointer :: p_KdiagonalA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A12 or NULL
    ! if not present
    real(DP), dimension(:), pointer :: p_DA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A21 or NULL
    ! if not present
    real(DP), dimension(:), pointer :: p_DA21 => null()


    ! Pointer to the matrix entries of the velocity matrix A44
    real(DP), dimension(:), pointer :: p_DA44 => null()

    ! Pointer to the matrix entries of the velocity matrix A55
    real(DP), dimension(:), pointer :: p_DA55 => null()


    ! Pointer to the matrix entries of the pressure identity matrix A33
    ! (if it exists).
    real(DP), dimension(:), pointer :: p_DA33 => null()

    ! Pointer to diagonal entries of A33
    integer, dimension(:), pointer :: p_KdiagonalA33 => null()

    ! Pointer to the matrix entries of the pressure identity matrix A66
    ! (if it exists).
    real(DP), dimension(:), pointer :: p_DA66 => null()

    ! Pointer to diagonal entries of A66
    integer, dimension(:), pointer :: p_KdiagonalA66 => null()


    ! Pointer to the column structure of the matrix A45 and A54
    integer, dimension(:), pointer :: p_KcolA45 => null()

    ! Pointer to the row structure of the matrix A45 and A54
    integer, dimension(:), pointer :: p_KldA45 => null()

    ! Pointer to diagonal entries in the matrix A45 and A54
    integer, dimension(:), pointer :: p_KdiagonalA45 => null()

    ! Pointer to the matrix entries of the velocity matrix A45 or NULL
    ! if not present
    real(DP), dimension(:), pointer :: p_DA45 => null()

    ! Pointer to the matrix entries of the velocity matrix A54 or NULL
    ! if not present
    real(DP), dimension(:), pointer :: p_DA54 => null()


    ! Pointer to the column structure of the mass matrix
    integer, dimension(:), pointer :: p_KcolM => null()

    ! Pointer to the row structure of the mass matrix
    integer, dimension(:), pointer :: p_KldM => null()

    ! Pointer to diagonal entries in the mass matrix
    integer, dimension(:), pointer :: p_KdiagonalM => null()

    ! Pointer to the matrix entries of the mass matrix at position
    ! (1,4) and (2,5) in the primal system, or NULL if not present.
    ! Has the same structure as the mass matrix.
    real(DP), dimension(:), pointer :: p_DM14 => null()
    real(DP), dimension(:), pointer :: p_DM25 => null()

    ! Pointer to the matrix entries of the mass matrix at position
    ! (1,5) and (2,4) in the primal system, or NULL if not present.
    ! Has the same structure as the mass matrix.
    real(DP), dimension(:), pointer :: p_DM15 => null()
    real(DP), dimension(:), pointer :: p_DM24 => null()

    ! Pointer to the coupling system at position (4,1), or NULL if not present
    ! Has the same structure as the mass matrix.
    real(DP), dimension(:), pointer :: p_DR41 => null()

    ! Pointer to the coupling system at position (5,2), or NULL if not present
    ! Has the same structure as the mass matrix.
    real(DP), dimension(:), pointer :: p_DR52 => null()

    ! Pointer to the coupling system at position (4,2), or NULL if not present
    ! Has the same structure as the mass matrix.
    real(DP), dimension(:), pointer :: p_DR42 => null()

    ! Pointer to the coupling system at position (5,1), or NULL if not present
    ! Has the same structure as the mass matrix.
    real(DP), dimension(:), pointer :: p_DR51 => null()


    ! Pointer to the column structure of the B/D-matrices.
    integer, dimension(:), pointer :: p_KcolB => null()

    ! Pointer to the row structure of the B/D-matrices
    integer, dimension(:), pointer :: p_KldB => null()

    ! Pointer to the entries of the B1-matrix
    real(DP), dimension(:), pointer :: p_DB1 => null()

    ! Pointer to the entries of the B2-matrix
    real(DP), dimension(:), pointer :: p_DB2 => null()

    ! Pointer to the entries of the D1-matrix
    real(DP), dimension(:), pointer :: p_DD1 => null()

    ! Pointer to the entries of the D2-matrix
    real(DP), dimension(:), pointer :: p_DD2 => null()


    ! Spatial discretisation structure for X-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrU => null()

    ! Spatial discretisation structure for Y-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrV => null()

    ! Spatial discretisation structure for pressure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrP => null()

    ! Multiplication factors for the submatrices; taken from the system matrix.
    ! (-> Not used in the current implementation! Although it is easy to include
    ! that into VANKA, some further speed analysis has to be done to make
    ! sure there is not too much speed impact when using these!)
    real(DP), dimension(6,6) :: Dmultipliers

  end type

!</typeblock>

!<typeblock>

  ! A structure that saves matrix pointers for the 3D-Navier-Stokes
  ! VANKA method for Navier-Stokes systems.

  type t_vankaPointer3DNavSt
    ! Pointer to the column structure of the velocity matrix A11, A22 and A33
    ! A11, A22 and A33 must have the same structure.
    integer, dimension(:), pointer :: p_KcolA => null()

    ! Pointer to the row structure of the velocity matrix A11, A22 and A33.
    ! They must have the same structure.
    integer, dimension(:), pointer :: p_KldA => null()

    ! Pointer to diagonal entries in the velocity matrix A11, A22 and A33
    integer, dimension(:), pointer :: p_KdiagonalA => null()

    ! Pointer to the matrix entries of the velocity matrix A11
    real(DP), dimension(:), pointer :: p_DA => null()

    ! Pointer to the matrix entries of the velocity matrix A22 or NULL
    ! if A11=A22.
    real(DP), dimension(:), pointer :: p_DA22 => null()

    ! Pointer to the matrix entries of the velocity matrix A33 or NULL
    ! if A11=A33.
    real(DP), dimension(:), pointer :: p_DA33 => null()

!    ! Pointer to the column structure of the velocity matrix A12 and A21
!    INTEGER, DIMENSION(:), POINTER :: p_KcolA12 => NULL()
!
!    ! Pointer to the row structure of the velocity matrix A12 and A21
!    INTEGER, DIMENSION(:), POINTER :: p_KldA12 => NULL()
!
!    ! Pointer to diagonal entries in the velocity matrix A12 and A21
!    INTEGER, DIMENSION(:), POINTER :: p_KdiagonalA12 => NULL()
!
!    ! Pointer to the matrix entries of the velocity matrix A12 or NULL
!    ! if not present
!    REAL(DP), DIMENSION(:), POINTER :: p_DA12 => NULL()
!
!    ! Pointer to the matrix entries of the velocity matrix A21 or NULL
!    ! if not present
!    REAL(DP), DIMENSION(:), POINTER :: p_DA21 => NULL()

    ! Pointer to the column structure of the B/D-matrices.
    integer, dimension(:), pointer :: p_KcolB => null()

    ! Pointer to the row structure of the B/D-matrices
    integer, dimension(:), pointer :: p_KldB => null()

    ! Pointer to the entries of the B1-matrix
    real(DP), dimension(:), pointer :: p_DB1 => null()

    ! Pointer to the entries of the B2-matrix
    real(DP), dimension(:), pointer :: p_DB2 => null()

    ! Pointer to the entries of the B3-matrix
    real(DP), dimension(:), pointer :: p_DB3 => null()

    ! Pointer to the entries of the D1-matrix
    real(DP), dimension(:), pointer :: p_DD1 => null()

    ! Pointer to the entries of the D2-matrix
    real(DP), dimension(:), pointer :: p_DD2 => null()

    ! Pointer to the entries of the D3-matrix
    real(DP), dimension(:), pointer :: p_DD3 => null()

    ! Pointer to the matrix entries of the pressure identity matrix A44
    ! (if it exists).
    real(DP), dimension(:), pointer :: p_DA44 => null()

    ! Pointer to diagonal entries of A33
    integer, dimension(:), pointer :: p_KdiagonalA44 => null()

    ! Spatial discretisation structure for X-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrU => null()

    ! Spatial discretisation structure for Y-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrV => null()

    ! Spatial discretisation structure for Z-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrW => null()

    ! Spatial discretisation structure for pressure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrP => null()

    ! Multiplication factors for the submatrices; taken from the system matrix.
    ! (-> Not used in the current implementation! Although it is easy to include
    ! that into VANKA, some further speed analysis has to be done to make
    ! sure there is not too much speed impact when using these!)
    real(DP), dimension(4,4) :: Dmultipliers

  end type

!</typeblock>

!<typeblock>

  ! A structure that represents a general VANKA configuration.

  type t_vanka

    ! The problem class, the VANKA should be applient to. One of the
    ! VANKAPC_xxxx constants, e.g. VANKAPC_2DNAVIERSTOKES.
    ! This constants decides where in the t_vankaXXXX structures below
    ! the actual data about the VANKA can be found.
    integer :: cproblemClass

    ! The subtype of VANKA that should handle the above problem class.
    ! One of the VANKATP_xxxx constants, e.g. VANKATP_DIAGONAL.
    integer :: csubtype

    ! This flag is set in the init-routine of VANKA and indicates whether
    ! an extended version of VANKA must be called that supports scaled
    ! matrices or a diagonal matrix in the pressure block.
    ! =0: standard type. =1: extended version
    integer :: csubsubtype = 0

    ! Configuration block with parameters for the general VANKA;
    ! only vaid if if cproblemClassVanka==VANKAPC_GENERAL.
    type(t_vankaGeneral) :: rvankaGeneral

    ! Configuration block with parameters for the 2D Navier-Stokes VANKA;
    ! only vaid if if cproblemClassVanka==VANKAPC_2DNAVIERSTOKES.
    type(t_vankaPointer2DNavSt) :: rvanka2DNavSt

    ! Configuration block with parameters for the 2D Navier-Stokes VANKA
    ! for optimal control problems;
    ! only vaid if if cproblemClassVanka==VANKAPC_2DNAVIERSTOKESOPTC.
    type(t_vankaPointer2DNavStOptC) :: rvanka2DNavStOptC

    ! Configuration block with parameters for the 3D Navier-Stokes VANKA;
    ! only vaid if if cproblemClassVanka==VANKAPC_3DNAVIERSTOKES.
    type(t_vankaPointer3DNavSt) :: rvanka3DNavSt

    ! -------------------- NEW IMPLEMENTATION --------------------

    ! Configuration block with parameters for the 2D Navier-Stokes VANKA;
    ! only vaid if if cproblemClassVanka==VANKAPC_NAVIERSTOKES2D.
    type(t_vanka_NavSt2D) :: rvankaNavSt2D

    ! Configuration block with parameters for the 2D Boussinesq VANKA;
    ! only vaid if if cproblemClassVanka==VANKAPC_BOUSSINESQ2D.
    type(t_vankaPointerBouss2D) :: rvankaBouss2D

    ! Configuration block with parameters for the 2D optimal control VANKA
    ! for Navier-Stokes; only vaid if if cproblemClassVanka==VANKAPC_NAVIERSTOKESOPTC2D.
    type(t_vanka_NavStOptC2D) :: rvankaNavStOptC2D

  end type

!</typeblock>

!</types>

!<constants>
!<constantblock description="Constants defining the blocking of element sets in VANKA">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of elements to handle simultaneously in general VANKA
#ifndef VANKA_NELEMSIM
  integer, parameter, public :: VANKA_NELEMSIM   = 1000
#endif

!</constantblock>
!</constants>



  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: vanka_perfconfig

  !************************************************************************

contains

  !****************************************************************************

!<subroutine>

  subroutine vanka_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      vanka_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(vanka_perfconfig)
      vanka_perfconfig%NELEMSIM = VANKA_NELEMSIM
    end if

  end subroutine vanka_initPerfConfig

  ! ***************************************************************************
  ! General VANKA for conformal discretisations.
  !
  ! The general VANKA is configured in vanka_initConformal to a special type
  ! problem (2D Navier Stokes e.g.) and applies the most suitable VANKA
  ! algorithm to a vector when being executed. Roughtly said, it is a wrapper
  ! for all types of VANKA that are realised in this module.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_initConformal (rmatrix,rvanka,cproblemClass,csubtype,rperfconfig)

!<description>
  ! Initialises the VANKA for conformal discretisations.
  ! Checks if the VANKA variant for conformal discretisations
  ! can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! VANKA does in general not support scaled matrices in rmatrix. The only
  ! exception is that by setting dscaleFactor=0.0 in one of the submatrices,
  ! a matrix can be deactivated.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in), target :: rmatrix

  ! The problem class that should be handled with this VANKA. One of the
  ! VANKAPC_xxxx constants, e.g. VANKAPC_2DNAVIERSTOKES
  integer, intent(in) :: cproblemClass

  ! The VANKA solver subtype that should handle the above problem class.
  ! One of the VANKATP_xxxx constants, e.g. VANKATP_DIAGONAL.
  integer, intent(in) :: csubtype

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vanka), intent(out) :: rvanka
!</output>

!</subroutine>

    select case (cproblemClass)
    case (VANKAPC_GENERAL)
      ! General VANKA
      call vanka_initGeneralVanka (rmatrix,rvanka%rvankaGeneral,rperfconfig)

    case (VANKAPC_2DNAVIERSTOKES)
      ! Vanka for 2D Navier-Stokes problems
      call vanka_init2DNavierStokes (rmatrix,rvanka,csubtype)

    case (VANKAPC_2DNAVIERSTOKESOPTC)
      ! Vanka for 2D Navier-Stokes optimal control problems
      call vanka_init2DNavierStokesOptC (rmatrix,rvanka)

    case (VANKAPC_3DNAVIERSTOKES)
      ! Vanka for 3D Navier-Stokes problems
      call vanka_init3DNavierStokes (rmatrix,rvanka)

    ! -------------------- NEW IMPLEMENTATION --------------------

    case (VANKAPC_NAVIERSTOKES2D)
      ! Vanka for 2D Navier-Stokes problems
      call vanka_init_NavSt2D(rvanka%rvankaNavSt2D, rmatrix, csubtype)

    case (VANKAPC_BOUSSINESQ2D)
      ! Vanka for 2D Boussinesq problems
      call vanka_initBoussinesq2D (rmatrix,rvanka%rvankaBouss2D,csubtype)

    case (VANKAPC_NAVIERSTOKESOPTC2D)
      ! Vanka for 2D Navier-Stokes optimal control problems
      call vanka_init_NavStOptC2D(rvanka%rvankaNavStOptC2D, rmatrix, csubtype)

    end select

    ! Initialise general data in the VANKA structure.
    rvanka%cproblemClass = cproblemClass
    rvanka%csubtype = csubtype

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_doneConformal (rvanka)

!<description>
  ! This routine cleans up a general VANKA solver. All memory allocated in
  ! rvankaGeneral is released.
!</description>

!<inputoutput>
  ! The VANKA structure to be cleaned up.
  type(t_vanka), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

    ! Type of VANKA?
    select case (rvanka%cproblemClass)
    case (VANKAPC_2DNAVIERSTOKES)
      ! Vanka for 2D Navier-Stokes problems
      call vanka_done2DNavierStokes (rvanka%rvanka2DNavSt)
    case (VANKAPC_GENERAL)
      ! Release data of the general VANKA
      call vanka_doneGeneralVanka (rvanka%rvankaGeneral)
    case (VANKAPC_NAVIERSTOKESOPTC2D)
      ! Vanka for 2D Navier-Stokes problems
      call vanka_done_NavStOptC2D (rvanka%rvankaNavStOptC2D)
    end select

    rvanka%cproblemClass = VANKAPC_NONE

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_conformal (rvanka, rvector, rrhs, domega, rperfconfig)

!<description>
  ! This routine applies the VANKA algorithm to the system $Ax=b$.
  ! This VANKA variant is rather general and can be applied to arbitrary
  ! conformal discretisations. It automatically chooses the correct
  ! VANKA 'subvariant' depending on the discretisation.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANKA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vanka), intent(inout) :: rvanka

  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! Which type of VANKA should be used?
    select case (rvanka%cproblemClass)
    case (VANKAPC_GENERAL)
      ! Use the general VANKA. This one handles the element distributions
      ! internally, so afterwards we can jump out of the loop without handling
      ! the other element distributions.

      call vanka_general (rvanka%rvankaGeneral, rvector, rrhs, domega, rperfconfig)

    case (VANKAPC_2DNAVIERSTOKES)
      ! 2D Navier Stokes problem.
      call vanka_2DNavierStokes (rvanka%rvanka2DNavSt, rvector, rrhs, domega,&
          rvanka%csubtype,rvanka%csubsubtype)

    case (VANKAPC_2DNAVIERSTOKESOPTC)
      ! 2D Navier Stokes problem.
      call vanka_2DNavierStokesOptC (rvanka%rvanka2DNavStOptC, rvector, rrhs, domega,&
          rvanka%csubtype)

    case (VANKAPC_3DNAVIERSTOKES)
      ! 3D Navier Stokes problem.
      call vanka_3DNavierStokes (rvanka%rvanka3DNavSt, rvector, rrhs, domega,&
          rvanka%csubtype,rvanka%csubsubtype)

    ! ---------------------- NEW IMPLEMENTATION ----------------------

    case (VANKAPC_NAVIERSTOKES2D)
      ! 2D Navier-Stokes problem.
      call vanka_solve_NavSt2D(rvanka%rvankaNavSt2D, rvector, rrhs, 1, domega)

    case (VANKAPC_BOUSSINESQ2D)
      ! 2D Boussinesq problem.
      call vanka_Boussinesq2D (rvanka%rvankaBouss2D, rvector, rrhs, domega,&
          rvanka%csubtype)

    case (VANKAPC_NAVIERSTOKESOPTC2D)
      ! 2D Navier-Stokes problem.
      call vanka_solve_NavStOptC2D(rvanka%rvankaNavStOptC2D, rvector, rrhs, 1, domega)

    case default
      call output_line ('Unknown VANKA problem class!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_conformal')
      call sys_halt()

    end select

  end subroutine

  ! ***************************************************************************
  ! GENERAL VANKA! Supports (non-transposed) matrices of any kind and size.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_initGeneralVanka (rmatrix,rvankaGeneral,rperfconfig)

!<description>
  ! This routine initialises the general VANKA solver and allocates
  ! necessary memory for the iteration.
!</description>

!<input>
  ! The system matrix that is used during the VANKA iteration.
  ! Remark: VANKA saves a pointer to this matrix, so the structure must exist
  !  until the system is solved! (Usually this points to the system matrix in
  !  the corresponding solver structure of the underlying linear solver...)
  type(t_matrixBlock), intent(in), target :: rmatrix

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! VANKA spiecific structure. Contains internal data and allocated memory.
  type(t_vankaGeneral), intent(out) :: rvankaGeneral
!</output>

!</subroutine>

    ! local variables
    logical :: bfirst
    integer :: nblocks,i,j,nmaxLocalDOFs,ndofsPerElement
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => vanka_perfconfig
    end if

    ! Make sure the block matrix is not rectangular!
    if (rmatrix%nblocksPerCol .ne. rmatrix%nblocksPerRow) then
      call output_line ('System matrix is rectangular!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initGeneralVanka')
      call sys_halt()
    end if

    nblocks = rmatrix%nblocksPerCol
    nmaxLocalDOFs = 0
    ndofsPerElement = 0

    ! Allocate memory for the matrix structures
    rvankaGeneral%nblocks = nblocks
    allocate(rvankaGeneral%p_Rmatrices(nblocks,nblocks))

    allocate(rvankaGeneral%p_InDofsLocal(nblocks))
    allocate(rvankaGeneral%p_IblockOffset(nblocks+1))

    ! This type of VANKA only supports a uniform discretisation
    ! and matrix format 7 or 9. Transposed matrices are not allowed.
    !
    ! Check all matrices that this is the case.
    ! Build the Rmatrices structure array. We manually create a
    ! (nblock,nblock) array inside of a 1-dimensional array and cast a rank
    ! change of the array on call to the actual VANKA subroutine later.
    !
    ! Loop through the columns of the block matrix.
    !
    ! Offset position of the first block is = 0.
    rvankaGeneral%p_IblockOffset(1) = 0

    do i=1,nblocks

      ! Note this block as 'not processed'
      bfirst = .true.

      ! Loop through the rows of the current matrix column.
      do j=1,nblocks
        if (lsysbl_isSubmatrixPresent(rmatrix,j,i)) then
          ! Get a/the discretisation structure of the current block/matrix column
          p_rdiscretisation => rmatrix%RmatrixBlock(j,i)%p_rspatialDiscrTrial

          if ((p_rdiscretisation%ccomplexity .ne. SPDISC_UNIFORM) .and. &
              (p_rdiscretisation%ccomplexity .ne. SPDISC_CONFORMAL)) then
            call output_line (&
                'General VANKA supports only uniform and conformal discretisations!',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_initGeneralVanka')
            call sys_halt()
          end if

          if ((rmatrix%RmatrixBlock(j,i)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
              (rmatrix%RmatrixBlock(j,i)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
            call output_line (&
                'General VANKA supports only matrix structure 7 and 9!',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_initGeneralVanka')
            call sys_halt()
          end if

          rvankaGeneral%p_Rmatrices(j,i)%bexists = .true.
          rvankaGeneral%p_Rmatrices(j,i)%dscaleFactor = &
            rmatrix%RmatrixBlock(j,i)%dscaleFactor
          rvankaGeneral%p_Rmatrices(j,i)%btransposed = &
            iand(rmatrix%RmatrixBlock(j,i)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0

          call lsyssc_getbase_double (rmatrix%RmatrixBlock(j,i),&
                                      rvankaGeneral%p_Rmatrices(j,i)%p_DA)
          call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(j,i),&
                                    rvankaGeneral%p_Rmatrices(j,i)%p_Kcol)
          call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(j,i),&
                                   rvankaGeneral%p_Rmatrices(j,i)%p_Kld)

          if (bfirst) then
            ! This block has not yet been processed.
            !
            ! Get the NEQ of the current block and save it as offset position
            ! in the global vector for the next block.
            rvankaGeneral%p_IblockOffset(i+1) = &
              rvankaGeneral%p_IblockOffset(i) + rmatrix%RmatrixBlock(j,i)%NCOLS

            ! We need some information for calculating DOF`s later.
            ! Get the number of local DOF`s in the current block.
            ! Note that we restrict to uniform discretisations!
            rvankaGeneral%p_InDofsLocal(i) = elem_igetNDofLoc(p_rdiscretisation% &
                                                RelementDistr(1)%celement)

            ! Calculate the maximum number of local DOF`s
            nmaxLocalDOFs = max(nmaxLocalDOFs,rvankaGeneral%p_InDofsLocal(i))

            ! Calculate the total number of local DOF`s
            ndofsPerElement = ndofsPerElement + rvankaGeneral%p_InDofsLocal(i)

            bfirst = .false.

          end if

        else
          rvankaGeneral%p_Rmatrices(j,i)%bexists = .false.
        end if

      end do
    end do

    if ((nmaxLocalDOFs .eq. 0) .or. (ndofsPerElement .eq. 0)) then
      call output_line (&
          'Matrix is empty!',OU_CLASS_ERROR,OU_MODE_STD,'vanka_initGeneralVanka')
      call sys_halt()
    end if

    ! Save the max. and total number of local DOF`s
    rvankaGeneral%nmaxLocalDOFs = nmaxLocalDOFs
    rvankaGeneral%ndofsPerElement = ndofsPerElement

    ! We know the maximum number of DOF`s now. For the later loop over the
    ! elements, allocate memory for storing the DOF`s of an element set.
    allocate(rvankaGeneral%p_IelementDOFs(nmaxLocalDOFs,p_rperfconfig%NELEMSIM,nblocks))

    ! Remember the matrix
    rvankaGeneral%p_rmatrix => rmatrix

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_doneGeneralVanka (rvankaGeneral)

!<description>
  ! This routine cleans up a general VANKA solver. All memory allocated in
  ! rvankaGeneral is released.
!</description>

!<inputoutput>
  ! The general-VANKA structure to be cleaned up.
  type(t_vankaGeneral), intent(inout) :: rvankaGeneral
!</inputoutput>

!</subroutine>

    ! Release memory allocated in the init-routine

    if (associated(rvankaGeneral%p_IelementDOFs)) &
      deallocate(rvankaGeneral%p_IelementDOFs)

    if (associated(rvankaGeneral%p_InDofsLocal)) &
      deallocate(rvankaGeneral%p_InDofsLocal)

    if (associated(rvankaGeneral%p_IblockOffset)) &
      deallocate(rvankaGeneral%p_IblockOffset)

    if (associated(rvankaGeneral%p_Rmatrices)) &
      deallocate(rvankaGeneral%p_Rmatrices)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_general (rvankaGeneral, rvector, rrhs, domega, rperfconfig)

!<description>
  ! This routine applies the general-VANKA algorithm to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvankaGeneral structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix A.
!</description>

!<input>

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The general-VANKA structure. Must have been initialised with
  ! vanka_initGeneralVanka before.
  type(t_vankaGeneral), intent(inout) :: rvankaGeneral

  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j
  integer :: IELmax, IELset, iel, ieldistr
  integer, dimension(:), pointer :: p_IelementList
  real(DP), dimension(:), pointer :: p_Drhs,p_Dvector
  integer, dimension(:), pointer :: p_Ipermutation
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

  if (present(rperfconfig)) then
    p_rperfconfig => rperfconfig
  else
    p_rperfconfig => vanka_perfconfig
  end if

  ! Saved matrix and the vector(s) must be compatible!
  call lsysbl_isMatrixCompatible(rvector,rvankaGeneral%p_rmatrix,.false.)
  call lsysbl_isVectorCompatible(rvector,rrhs)

  ! Get the data arrays of the vector/rhs
  call lsysbl_getbase_double (rvector,p_Dvector)
  call lsysbl_getbase_double (rrhs,p_Drhs)

  ! Get the discretisation structure that tells us which elements form
  ! element groups...
  p_rdiscretisation => rvector%RvectorBlock(1)%p_rspatialDiscr

  ! Loop over the element distributions/groups
  do ieldistr = 1,p_rdiscretisation%inumFEspaces

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions.
    call storage_getbase_int (p_rdiscretisation% &
                              RelementDistr(ieldistr)%h_IelementList, &
                              p_IelementList)

    ! Loop over the elements - blockwise.
    do IELset = 1, size(p_IelementList), p_rperfconfig%NELEMSIM

      ! We always handle NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most NELEMSIM
      ! elements simultaneously.
      IELmax = min(size(p_IelementList),IELset-1+p_rperfconfig%NELEMSIM)

      ! Loop over the blocks in the block vector to get the DOF`s everywhere.

      do i=1,rvector%nblocks

        ! Calculate the global DOF`s of all blocks.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rvector%RvectorBlock(i)%p_rspatialDiscr,&
                                     p_IelementList(IELset:IELmax), &
                                     rvankaGeneral%p_IelementDOFs(:,:,i))

        ! If the vector is sorted, push the DOF`s through the permutation to get
        ! the actual DOF`s.
        if (rvector%RvectorBlock(i)%isortStrategy .gt. 0) then

          call storage_getbase_int(rvector%RvectorBlock(i)%h_IsortPermutation,&
                                  p_Ipermutation)

          do iel=1,IELmax-IELset+1
            do j=1,rvankaGeneral%p_InDofsLocal(i)
              ! We are not resorting the vector but determining the 'sorted'
              ! DOF`s - this needs the 2nd half of the permutation.
              rvankaGeneral%p_IelementDOFs(j,iel,i) = &
                 p_Ipermutation(rvankaGeneral%p_IelementDOFs(j,iel,i) &
                +rvector%RvectorBlock(i)%NEQ)
            end do
          end do
        end if

      end do

      ! Now, IdofsTotal contains all DOF`s on each element, over all discretisations.
      !
      ! Call the actual VANKA to process the DOF`s on each element.
      call vanka_general_double_mat79 (p_Dvector, p_Drhs, domega, &
          rvankaGeneral%p_Rmatrices,IELmax-IELset+1,&
          rvankaGeneral%p_IblockOffset,rvankaGeneral%nblocks,&
          rvankaGeneral%p_InDofsLocal,rvankaGeneral%ndofsPerElement,&
          rvankaGeneral%p_IelementDOFs)

    end do

  end do ! ieldistr

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_general_double_mat79 (Dvector, Drhs, domega, Rmatrices,&
             nelements,IblockOffset,nblocks,InDofsLocal,ndofsPerElement,&
             IelementDofs)

!<description>
  ! This routine applies one step of the VANKA solver (local block
  ! gauss-seidel) to a given set of solution/RHS vector.
  ! All given matrices are assumed to be format 7/9,
  ! double precision. The given vector is assumed to be double precision.
!</description>

!<input>
  ! The (block) RHS vector, given as one large array.
  real(DP), dimension(:), intent(in) :: Drhs

  ! A relaxation parameter. Standard = 1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of matrices to handle; directly specified by pointers
  ! to the substructures (data/columns/rows).
  type(t_matrixPointer79Vanka), dimension(:,:),&
                                intent(in) :: Rmatrices

  ! Number of elements that should be processed in this sweep.
  integer, intent(in) :: nelements

  ! Number of blocks in the vectors
  integer, intent(in) :: nblocks

  ! Offset position of the blocks in the vector.
  ! Block i starts at position IblockOffset(i)+1 in Dvector / Drhs.
  ! IblockOffset(nblocks+1) gives the number of equations in Dvector/Drhs.
  integer, dimension(:), intent(in) :: IblockOffset

  ! Number of local DOF`s in each block.
  integer, dimension(:), intent(in) :: InDofsLocal

  ! Total number of local DOF`s per element
  integer, intent(in) :: ndofsPerElement

  ! List of DOF`s on every element for every block.
  ! DIMENSION(nmaxDOFs,nelements,nblocks)
  integer, dimension(:,:,:), intent(in) :: IelementDOFs

!</input>

!<inputoutput>
  ! The initial (block) solution vector. Is overwritten by the new (block)
  ! solution vector.
  real(DP), dimension(:), intent(inout) :: Dvector
!</inputoutput>

!</subroutine>

    ! One iteration of vanka smother (block gauss-seidel) on the system
    !
    !    A11 A12 ... A1n | U1    F1
    !    A21 A22 ... A2n | U2  = F2
    !     :   :  ...  :  |  :     :
    !    An1 An2 ... Ann | Un  = Fn
    !
    ! Let us first describe the method used here.
    !
    ! The above block system can be written as a general block system of
    ! the form
    !
    !             A x = f
    !
    ! Consider a general defect correction approach for this system:
    !
    !     x_{n+1}  =  x_n  +  \omega C^{-1} (b - A x_n)
    !
    ! With some damping parameter \omega and some preconditioner C^{-1}.
    ! Normally, this iteration is used globally - and the usual linear
    ! solver that calls this routine usually realises this approach.
    ! VANKA now applies this defect correction loop *locally*, i.e.
    ! not for the full system at once.
    !
    ! This local approach is based on a geometric point of view.
    ! In general, one could imagine a set of connected cells in the
    ! global domain \Omega where to apply the algorithm to (the LMPSC
    ! approach), but for simplicity, consider only one cell. Again for
    ! simplicity, imagine that our FE-spaces is Q1~/Q0.
    !
    ! We loop over each cell in the domain, one after the other, and
    ! change the DOF`s in the solution vector there:
    !
    ! We fetch all the data (e.g. velocity, pressure) on that cell. On the
    ! first cell, we have only "old" velocity entries. These values
    ! are updated and then the calculation proceeds with the 2nd cell.
    !
    !        old                      new
    !     |---X---|                |---X---|
    !     |       |                |       |
    ! old X   1   X old   -->  new X   1   X new
    !     |       |                |       |
    !     |---X---|                |---X---|
    !        old                      new
    !
    ! From the second cell on, there might be "old" data and "new"
    ! data on that cell - the old data that has not been updated and
    ! perhaps some already updated velocity data from a neighbor cell.
    !
    !        new     old                   new     new
    !     |---X---|---X---|             |---X---|---X---|
    !     |       |       |             |       |       |
    ! new X   1   X   2   X old --> new X   1   X   2   X new
    !     |       |new    |             |       |newer  |
    !     |---X---|---X---|             |---X---|---X---|
    !        new     old                   new     new
    !
    ! These values are updated and then the calculation proceeds
    ! with the next cell.
    ! As can be seen in the above picture, the "new" node in the
    ! middle is even going to be a "newer" node when handled again
    ! for the 2nd cell. This is meant by "Gauss-Seldel" character:
    ! Information is updated subsequently by using "old" data and
    ! "new" data from a previous calculation.
    !
    ! So how to do this in detail?
    !
    ! Again consider the problem:
    !
    !             A x = f
    !
    ! We assume, that all components in the vector x are
    ! given - except for the unknowns current element; these unknowns
    ! are located anywhere in the vector x. The idea is to
    ! shift "everything known" to the right hand side to obtain
    ! a system for only these unknowns!
    !
    ! We extract all the lines of the system that correspond to
    ! our 'unknown' DOF`s on our element; this results in a rectangular
    ! system of the form
    !
    !    [ === A~ === ] x  = (f~)
    !
    ! So #rows(A~)=#rows(f~)=#DOF`s on the element! Furthermore we set
    ! Now we make a defect-correction approach for this system:
    !
    !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
    !                                     -----------
    !                                        =d~
    !
    ! Here the 'projection' operator simply converts the small
    ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
    ! of the same size as x - what is easy using the number of
    ! the DOF`s on the element.
    !
    ! The only question now will be: What is C^{-1}?
    !
    ! Well, here we have different choices. Everything depends on the
    ! matrix A~, which is unfortunately rectangular: It is a
    ! (#DOF`s on the element, #DOF`s in the space) matrix - so
    ! in case of a Q1~/Q0 discretisation, it is a (5,NEQ) matrix!
    !
    ! For full linear systems, one would choose C=A, which ist the
    ! theoretically best preconditioner. What we do here is simply
    ! extracting all columns of A~ that correspond to the DOF`s
    ! on the current element:
    !
    !   C:=delete columns of A~ that do not belong to DOF`s on the element
    !
    ! This then leads to a square preconditioner C^{-1} - and that is the
    ! full method, because C^{-1} can be applied directly using Lapack e.g.!
    !
    !
    ! So all in all we write our algorithm in short:
    !
    ! loop over all elements in the given element set
    !   extract local matrix aa
    !   extract the local rhs ff
    !   compute local residuum: ff := ff - aa * x
    !   solve: xx := aa^(-1) ff
    !   update global solution vector: x := x + omega*xx
    ! end loop

    ! local variables
    integer :: iel
    integer, dimension(max(nblocks,1)+1) :: IlocalIndex
    integer, dimension(ndofsPerElement) :: IlocalDOF,IglobalDOF
    integer :: i,j,k,iidx,iminiDOF
    integer :: irow,idof
    integer :: icol

    ! Memory for our local system; let us hope it is not too big :)
    real(DP), dimension(ndofsPerElement,ndofsPerElement) :: Daa
    real(DP), dimension(ndofsPerElement) :: Dff
    real(DP) :: dscale
    integer, dimension(ndofsPerElement) :: Ipiv
    integer :: iinfo

    ! Quickly check the matrices if one of them is saved transposed.
    ! VANKA does not support transposed matrices; would kill computational
    ! time, as we cannot extract columns from a structure-9 matrix with
    ! reasonable effort!
    do i=1,size(Rmatrices,1)
      do j=1,size(Rmatrices,2)
        if (Rmatrices(i,j)%bexists .and. Rmatrices(i,j)%btransposed) then
          call output_line (&
              'General VANKA does not support transposed matrices!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_general_double_mat79')
          call sys_halt()
        end if
      end do ! j
    end do ! i

    ! Build an index pointer for accessing local DOF`s
    IlocalIndex(1) = 0
    do i=2,nblocks+1
      IlocalIndex(i) = IlocalIndex(i-1)+InDOFsLocal(i-1)
    end do

    ! Ok, let us start with the loop over the elements in the given
    ! element list.
    do iel = 1,nelements

      ! IelementDOFs (.,iel,.) gives now for every block in the system
      ! the DOF`s on this element.
      !
      ! First copy the RHS entries of f to f~.

      iidx = 1
      do i=1,nblocks
        do j=1,InDOFsLocal(i)
          iidx = IlocalIndex(i)+j

          ! Get the DOF on the element:
          idof = IelementDOFs(j,iel,i)

          ! Number of DOF relative to this block
          IlocalDOF(iidx) = idof

          ! Calculate corresponding global DOF:
          IglobalDOF(iidx) = idof+IblockOffset(i)

          Dff(iidx) = Drhs(IglobalDOF(iidx))
        end do
      end do

      ! Compute  ff := ff - A x  for the local unknowns
      ! to build the local residuum Dff = f~-A~x.
      !
      ! Simultaneously extract the local matrix into one array Daa(:,:).
      ! But first clear the matrix - maybe that there are unused subblocks!
      Daa = 0

      ! Loop through the rows of the block matrix
      do i=1,size(Rmatrices,1)

        ! Loop through the columns of the block matrix
        do j=1,size(Rmatrices,2)

          ! Is there a matrix saved at this position?
          if (Rmatrices(i,j)%bexists .and. &
              Rmatrices(i,j)%dscaleFactor .ne. 0.0_DP) then

            dscale = Rmatrices(i,j)%dscaleFactor

            ! Block column j in the global matrix corresponds to
            ! block column j in the small local matrix we have to fill with data.
            !
            !   ....      :
            !   ....      :           a11 a12         a15
            !                         a21 a22         a25
            !        .... :   ==>             a33 a34 a35
            !        .... :                   a43 a44 a45
            !   .... ....             a51 a52 a53 a54
            !   ==== ==== =           ======= ======= ===
            !    1    2   3 corr. to     1       2     3
            !
            !     n x n-mat.          5x5-mat or so
            !
            ! From IlocalIndex, get the starting address of the j-th block
            ! in the local solution vector. In the above example, for j=2 e.g.
            ! this gives the starting address of the two global DOF`s
            ! that correspond to the columns 3 and 4 in the local matrix.
            iidx = IlocalIndex(j)

            ! Loop through the DOF`s that correspond to this block:
            do irow = IlocalIndex(i)+1,IlocalIndex(i+1)

              ! Get the actual DOF, relative to this block.
              ! This is the row of the matrix that must be multiplied by x
              ! and be subtracted from the local RHS.
              idof = IlocalDOF(irow)

              ! Loop through the row of the matrix to its contribution to "b-Ax".
              do k = Rmatrices(i,j)%p_Kld(idof) , Rmatrices(i,j)%p_Kld(idof+1)-1

                ! Get the column number in the global matrix. This gives
                ! the global DOF = number of the element in x that must be multiplied
                ! with that matrix entry.
                icol = Rmatrices(i,j)%p_Kcol(k)+IblockOffset(j)

                ! Build the defect
                Dff(irow) = Dff(irow) &
                          - dscale * Rmatrices(i,j)%p_Da(k) * Dvector(icol)

                ! icol is the number of a DOF.
                ! Check if this DOF belongs to the DOF`s we have to
                ! extract to our local system.
                ! Loop through the DOF`s corresponding to column j of the block system.
                ! In the above example, this checks only the two DOF`s corresponding
                ! to a?3 and a?4.
                do iminiDOF = 1,InDOFsLocal(j)
                  if (icol .eq. IglobalDOF(iidx+iminiDOF)) then
                    ! Yes. Get the matrix entry, write it to the local matrix
                    Daa(irow,iidx+iminiDOF) = dscale*Rmatrices(i,j)%p_Da(k)
                    exit
                  end if
                end do ! idof

              end do ! k

            end do ! irow

          end if ! exists

        end do ! j

      end do ! i

      ! Ok, we now have our local matrix and our local defect vector d~.
      ! Apply LAPACK to the local system to solve C^{-1} d~.

      call DGETRF( ndofsPerElement, ndofsPerElement, Daa, ndofsPerElement, &
                   Ipiv, iinfo )

      ! Note: It may happen that the matrix is singular!
      !  That is the case if all DOF`s are Dirichlet-DOF`s - for example
      !  if the element is completely inside of a rigid fictitious boundary
      !  object.
      ! What to do in this case? Nothing! Ignore the system!
      ! Why?
      !  - The values for the 'velocity' DOF`s are fixed, so it is no use
      !    to try to calculate them.
      !  - If there is a zero-block involved in a saddle-point problem,
      !    the pressure (that corresponds to the zero-block) is not
      !    connected to the velocity - it is completely free!
      !    By ignoring this system, we let the pressure as it is.
      !
      ! One can theoretically also solve a least-squares system with
      ! LAPACK. This calculate an y with ||f~ - C y|| -> min and |y|->min.
      ! But this is (because of the 0-block in A and therefore also in C
      ! and because of the unit-vectors in the other rows of C)
      ! exactly the case if the 'velocity' y matches f~ and the 'pressure'
      ! components are zero - which means nothing else than that there is
      ! no contribution in the preconditioned defect to correct the
      ! pressure entries of x.

      if (iinfo .eq. 0) then
        call DGETRS('N', ndofsPerElement, 1, Daa, ndofsPerElement, &
                    Ipiv, Dff, ndofsPerElement, iinfo )
        if (iinfo .eq. 0) then

          ! Dff is in-situ replaced by the solution - the preconditioned
          ! defect. Add this to our original solution vector to complete
          ! the 'local' defect-correction.

          do i=1,ndofsPerElement
            j = IglobalDOF(i)
            Dvector(j) = Dvector(j) + domega*Dff(i)
          end do

        end if
      end if

    end do ! iel

  end subroutine

! *****************************************************************************
! Problem class: VANKA variants for 2D Navier-Stokes problems
! *****************************************************************************

!<subroutine>

  subroutine vanka_init2DNavierStokes (rmatrix,rvanka,csubtype)

!<description>
  ! Initialises the VANKA variant for 2D Navier-Stokes problems
  ! for conformal discretisations.
  ! Checks if the "2D-Navier-Stokes" VANKA variant
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! The substructure rvanka%rvanka2DNavSt is intitialised according
  ! to the information provided in rmatrix.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in), target :: rmatrix

  ! Desired subtype
  integer, intent(in) :: csubtype
!</input>

!<inputoutput>
  ! t_vankaPointer2DSPNavSt structure that saves algorithm-specific parameters.
  type(t_vanka), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

    integer :: i,j
    logical :: bextended
    type(t_blockDiscretisation), pointer :: p_rblockDiscr

    bextended = .false.

    ! Matrix must be 3x3.
    if ((rmatrix%nblocksPerCol .ne. 3) .or. (rmatrix%nblocksPerRow .ne. 3)) then
      call output_line ('System matrix is not 3x3.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
      call sys_halt()
    end if

    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(3,:) must be (virtually) transposed. All matrices must be double precision.
    do i=1,3
      do j=1,3

        if (lsysbl_isSubmatrixPresent(rmatrix,i,j)) then

          if (i .le. 2) then
            if (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .ne. 0) then
              call output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
              call sys_halt()
            end if
          else
            if ((i .le. 2) .and. &
               (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .eq. 0)) then
              call output_line ('B1/B2 submatrices must be virtually',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
              call output_line ('transposed (LSYSSC_MSPEC_TRANSPOSED)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
              call sys_halt()
            end if
          end if

          if ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX9)) then
            call output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
            call sys_halt()
          end if

          if (rmatrix%RmatrixBlock(i,j)%cdataType .ne. ST_DOUBLE) then
            call output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
            call sys_halt()
          end if

          ! For scaled matrices, we have to use an extended sub-version of VANKA.
          if ((rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 1.0_DP) .and. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 0.0_DP)) then
            bextended = .true.
          end if

          if ((i .eq. j) .and. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 1.0_DP) ) then
            bextended = .true.
          end if

        end if ! neq != 0
      end do
    end do

    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,3)%NA .ne. rmatrix%RmatrixBlock(3,1)%NA) .or. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .ne. rmatrix%RmatrixBlock(3,1)%NCOLS)) then
      call output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,3)%NA .ne. rmatrix%RmatrixBlock(3,2)%NA) .or. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .ne. rmatrix%RmatrixBlock(3,2)%NCOLS)) then
      call output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
      call sys_halt()
    end if

    ! Fill the output structure with data of the matrices.
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanka%rvanka2DNavSt%p_DA )
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),&
        rvanka%rvanka2DNavSt%p_DB1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),&
        rvanka%rvanka2DNavSt%p_DB2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),&
        rvanka%rvanka2DNavSt%p_DD1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),&
        rvanka%rvanka2DNavSt%p_DD2)
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),&
        rvanka%rvanka2DNavSt%p_KcolB)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), &
        rvanka%rvanka2DNavSt%p_KldB )
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanka%rvanka2DNavSt%p_KcolA)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanka%rvanka2DNavSt%p_KldA )
    if (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanka%rvanka2DNavSt%p_KdiagonalA)
    else
      rvanka%rvanka2DNavSt%p_KdiagonalA => rvanka%rvanka2DNavSt%p_KldA
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,3,3)) then

      ! The matrix must be of format 7 or 9.
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
          rvanka%rvanka2DNavSt%p_DA33 )

      if (rmatrix%RmatrixBlock(3,3)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(3,3), &
                                rvanka%rvanka2DNavSt%p_KdiagonalA33)
      else
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,3), &
                                rvanka%rvanka2DNavSt%p_KdiagonalA33)
      end if

      ! The presence of A(3,3) forces the extended VANKA to be used
      bextended = .true.

    end if

    if (bextended) rvanka%csubsubtype = 1

    ! What is with A22? Is it the same as A11?
    if (.not. lsyssc_isMatrixContentShared (&
        rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2)) ) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
          rvanka%rvanka2DNavSt%p_DA22 )
    end if

    ! What is with A12 and A21? Do they exist? With a scale factor = 1.0?
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2) .and. &
        (rmatrix%RmatrixBlock(1,2)%dscaleFactor .eq. 1.0_DP)) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
          rvanka%rvanka2DNavSt%p_DA12 )

      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,2),&
          rvanka%rvanka2DNavSt%p_KcolA12)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,2), &
          rvanka%rvanka2DNavSt%p_KldA12 )

      ! Get the structure. It is assumed that A12 and A21 have the same!
      if (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,2), &
                                rvanka%rvanka2DNavSt%p_KdiagonalA12)
      else
        rvanka%rvanka2DNavSt%p_KdiagonalA12 => rvanka%rvanka2DNavSt%p_KldA12
      end if

      if (.not. lsysbl_isSubmatrixPresent(rmatrix,2,1)) then
        call output_line ('If A12 is given, A21 must also be given!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesO')
        call sys_halt()
      end if

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
          rvanka%rvanka2DNavSt%p_DA21 )
    end if

    ! Get the multiplication factors of the submatrices.
    rvanka%rvanka2DNavSt%Dmultipliers(1:3,1:3) = &
        rmatrix%RmatrixBlock(1:3,1:3)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest

    if (.not. associated(p_rblockDiscr)) then
      call output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
      call sys_halt()
    end if

    ! Get the discretisation structure of U,V and P from the block
    ! discretisation structure.
    rvanka%rvanka2DNavSt%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanka%rvanka2DNavSt%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanka%rvanka2DNavSt%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)

    if (rvanka%rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces .ne. &
        rvanka%rvanka2DNavSt%p_rspatialDiscrV%inumFESpaces) then
      call output_line (&
          'Discretisation structures of X- and Y-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
      call sys_halt()
    end if

    if ((rvanka%rvanka2DNavSt%p_rspatialDiscrP%inumFESpaces .ne. 1) .and. &
        (rvanka%rvanka2DNavSt%p_rspatialDiscrP%inumFESpaces .ne. &
          rvanka%rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces)) then
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      call output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokes')
      call sys_halt()
    end if

    ! Solution based VANKA variants need an additional temp vector
    if (csubtype .eq. VANKATP_DIAGONAL_SOLBASED) then
      call lsysbl_createVecBlockIndMat (rmatrix,rvanka%rvanka2DNavSt%rtempVector, .false.)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_done2DNavierStokes (rvanka2DNavSt)

!<description>
  ! Releases any memory allocated in rvanka2DNavSt.
!</description>

!<inputoutput>
  ! t_vanka structure that to be cleaned up.
  type(t_vankaPointer2DNavSt), intent(inout) :: rvanka2DNavSt
!</inputoutput>

!</subroutine>

    if (rvanka2DNavSt%rtempVector%NEQ .ne. 0) then
      call lsysbl_releaseVector (rvanka2DNavSt%rtempVector)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DNavierStokes (rvanka2DNavSt, rvector, rrhs, domega, &
      csubtype, csubsubtype)

!<description>
  ! This routine applies the VANKA variant for 2D Navier-Stokes problems
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANKA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! The subtype of VANKA that should handle the above problem class.
  ! One of the VANKATP_xxxx constants, e.g. VANKATP_DIAGONAL.
  integer :: csubtype

  ! The sub-subtype of VANKA that should handle the above problem class.
  ! =0: use standard VANKA. =1: use extended VANKA (e.g. with different
  !     multipliers in the matrices)
  integer :: csubsubtype

!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(inout) :: rvector

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(inout) :: rvanka2DNavSt
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ielementdist
    integer, dimension(:), pointer :: p_IelementList
    type(t_elementDistribution), pointer :: p_relementDistrU
    type(t_elementDistribution), pointer :: p_relementDistrV
    type(t_elementDistribution), pointer :: p_relementDistrP

    ! 2D Navier Stokes problem.

    ! Loop through the element distributions of the velocity.
    do ielementdist = 1,rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces

      ! Get the corresponding element distributions of U, V and P.
      p_relementDistrU => &
          rvanka2DNavSt%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanka2DNavSt%p_rspatialDiscrV%RelementDistr(ielementdist)

      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      if (rvanka2DNavSt%p_rspatialDiscrP%inumFESpaces .gt. 1) then
        p_relementDistrP => &
            rvanka2DNavSt%p_rspatialDiscrP%RelementDistr(ielementdist)
      else
        p_relementDistrP => &
            rvanka2DNavSt%p_rspatialDiscrP%RelementDistr(1)
      end if

      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y-velocity (and to that
      ! of the pressure).
      call storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)

      ! Which element combination do we have now?
      if ((elem_getPrimaryElement(p_relementDistrU%celement) .eq. EL_Q1T) .and. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .eq. EL_Q1T) .and. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .eq. EL_Q0)) then
        ! Q1~/Q1~/Q0 discretisation

        ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
        select case (csubtype)
        case (VANKATP_DIAGONAL)
          ! Diagonal VANKA; check if A12 exists.
          if (.not. associated(rvanka2DNavSt%p_DA12)) then
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            if (rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! Uniform discretisation
              call vanka_2DSPQ1TQ0simple (rvanka2DNavSt, &
                  rvector, rrhs, domega)
            else
              ! Conformal discretisation
              call vanka_2DSPQ1TQ0simpleConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          else
            ! Apply the conformal VANKA that allows different matrices
            ! in A11, A12, A21 and A22!
            call vanka_2DSPQ1TQ0simpleCoupConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
          end if

        case (VANKATP_DIAGONAL_SOLBASED)
          ! Diagonal VANKA. Solution based if possible. Check if A12 exists.
          if (.not. associated(rvanka2DNavSt%p_DA12)) then
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            if (rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! Uniform discretisation
              call vanka_2DSPQ1TQ0simpleSol (rvanka2DNavSt, &
                  rvector, rrhs, domega,rvanka2DNavSt%rtempVector)
            else
              ! Conformal discretisation
              call vanka_2DSPQ1TQ0simpleConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          else
            ! Apply the conformal VANKA that allows different matrices
            ! in A11, A12, A21 and A22!
            call vanka_2DSPQ1TQ0simpleCoupConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
          end if

        case (VANKATP_FULL)
          ! Full VANKA; check if A12 exists.
          if (.not. associated(rvanka2DNavSt%p_DA12)) then
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            ! Note: Up to now, there is no 'full' variant -- has to be implemented!
            if (rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! uniform discretisation;
              ! here, use the same as for the general conformal discretisation.
              ! Could be speeded up by introducing another variant...
              call vanka_2DSPQ1TQ0fullConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            else
              ! Conformal discretisation
              call vanka_2DSPQ1TQ0fullConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          else
            ! Apply the conformal VANKA that allows different matrices
            ! in A11, A12, A21 and A22!
            ! If we have multiplication factors, we even have to use an extended
            ! version of this.
            if (csubsubtype .eq. 0) then
              call vanka_2DSPQ1TQ0fullCoupConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            else
              call vanka_2DNSQ1TQ0fullCoupConfExt (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          end if

        case default
          call output_line (&
              'Unknown VANKA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNavierStokes')
          call sys_halt()

        end select

      ! Which element combination do we have now?
      else if ((elem_getPrimaryElement(p_relementDistrU%celement) .eq. EL_P1T) .and. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .eq. EL_P1T) .and. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .eq. EL_P0)) then
        ! Q1~/Q1~/Q0 discretisation

        ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
        select case (csubtype)
        case (VANKATP_DIAGONAL,VANKATP_DIAGONAL_SOLBASED)
          ! Diagonal VANKA; check if A12 exists.
          if (.not. associated(rvanka2DNavSt%p_DA12)) then
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            call vanka_2DSPP1TP0simpleConf (rvanka2DNavSt, &
                rvector, rrhs, domega,p_IelementList)
          else
            ! Apply the conformal VANKA that allows different matrices
            ! in A11, A12, A21 and A22!
            call vanka_2DSPP1TP0simpleCoupConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
          end if

        case (VANKATP_FULL)
          ! Full VANKA; check if A12 exists.
          if (.not. associated(rvanka2DNavSt%p_DA12)) then
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            ! Note: Up to now, there is no 'full' variant -- has to be implemented!
            if (rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! uniform discretisation;
              ! here, use the same as for the general conformal discretisation.
              ! Could be speeded up by introducing another variant...
              call vanka_2DSPP1TP0fullConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            else
              ! Conformal discretisation
              call vanka_2DSPP1TP0fullConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          else
            ! Apply the conformal VANKA that allows different matrices
            ! in A11, A12, A21 and A22!
            ! If we have multiplication factors, we even have to use an extended
            ! version of this.
            if (csubsubtype .eq. 0) then
              call vanka_2DSPP1TP0fullCoupConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            else
              call vanka_2DNSP1TP0fullCoupConfExt (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          end if

        case default
          call output_line (&
              'Unknown VANKA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNavierStokes')
          call sys_halt()

        end select

      else if &
        ((elem_getPrimaryElement(p_relementDistrU%celement) .eq. EL_Q2) .and.&
          (elem_getPrimaryElement(p_relementDistrV%celement) .eq. EL_Q2) .and.&
          (elem_getPrimaryElement(p_relementDistrP%celement) .eq. EL_QP1)) then
        ! Q2/Q2/QP1 discretisation

        ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
        select case (csubtype)
        case (VANKATP_DIAGONAL)
          ! Diagonal VANKA; check if A12 exists.
          if (.not. associated(rvanka2DNavSt%p_DA12)) then
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            if (rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! Uniform discretisation
              call vanka_2DSPQ2QP1simple (rvanka2DNavSt, &
                  rvector, rrhs, domega)
            else
              ! Conformal discretisation
              call vanka_2DSPQ2QP1simpleConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          else
            ! Apply the conformal VANKA that allows different matrices
            ! in A11, A12, A21 and A22!
            call vanka_2DSPQ2QP1simpleCoupConf (rvanka2DNavSt, &
                rvector, rrhs, domega,p_IelementList)
          end if

        case (VANKATP_FULL)
          ! Full VANKA; check if A12 exists.
          if (.not. associated(rvanka2DNavSt%p_DA12)) then
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            if (rvanka2DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! Uniform discretisation
              call vanka_2DSPQ2QP1full (rvanka2DNavSt, &
                  rvector, rrhs, domega)
            else
              ! Conformal discretisation
              call vanka_2DSPQ2QP1fullConf (rvanka2DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
          else
            ! Apply the conformal VANKA that allows different matrices
            ! in A11, A12, A21 and A22!
            call vanka_2DSPQ2QP1fullCoupConf (rvanka2DNavSt, &
                rvector, rrhs, domega,p_IelementList)
          end if

        case default
          call output_line (&
              'Unknown VANKA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNavierStokes')
          call sys_halt()

        end select

      else
        call output_line (&
            'Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNavierStokes')
        call sys_halt()

      end if

    end do

  end subroutine

  ! ***************************************************************************
  ! 2D Navier-Stokes VANKA, simple diagonal version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_init2DSPQ1TQ0simple (rmatrix,rvanka)

!<description>
  ! Checks if the "2D-Navier-Stokes-Q1T-Q0" VANKA variant can be applied to
  ! the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>

!<output>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(out) :: rvanka
!</output>

!</subroutine>

    integer :: i,j

    ! Matrix must be 3x3.
    if ((rmatrix%nblocksPerCol .ne. 3) .or. (rmatrix%nblocksPerRow .ne. 3)) then
      call output_line ('System matrix is not 3x3.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
      call sys_halt()
    end if

    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(3,:) must be (virtually) transposed. All matrices must be double precision.
    do i=1,3
      do j=1,3

        if (lsysbl_isSubmatrixPresent(rmatrix,i,j)) then

          if (i .le. 2) then
            if (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .ne. 0) then
              call output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
              call sys_halt()
            end if
          else
            if (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .eq. 0) then
              call output_line ('B1/B2 submatrices must be virtually '//&
                  'transposed (LSYSSC_MSPEC_TRANSPOSED)',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
              call sys_halt()
            end if
          end if

          if ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX9)) then
            call output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
            call sys_halt()
          end if

          if (rmatrix%RmatrixBlock(i,j)%cdataType .ne. ST_DOUBLE) then
            call output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
            call sys_halt()
          end if

          if (rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 1.0_DP) then
            call output_line ('Scaled matrices not supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
            call sys_halt()
          end if

        end if ! neq != 0
      end do
    end do

    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,3)%NA .ne. rmatrix%RmatrixBlock(3,1)%NA) .or. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .ne. rmatrix%RmatrixBlock(3,1)%NCOLS)) then
      call output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,3)%NA .ne. rmatrix%RmatrixBlock(3,2)%NA) .or. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .ne. rmatrix%RmatrixBlock(3,2)%NCOLS)) then
      call output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DSPQ1TQ0simple')
      call sys_halt()
    end if

    ! Fill the output structure with data of the matrices.
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),rvanka%p_DA )
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),rvanka%p_DB1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),rvanka%p_DB2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),rvanka%p_DD1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),rvanka%p_DD2)
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),rvanka%p_KcolB)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), rvanka%p_KldB )
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),rvanka%p_KcolA)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), rvanka%p_KldA )
    if (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                               rvanka%p_KdiagonalA)
    else
      rvanka%p_KdiagonalA => rvanka%p_KldA
    end if

    ! Get the multiplication factors of the submatrices
    rvanka%Dmultipliers(:,:) = rmatrix%RmatrixBlock(1:3,1:3)%dscaleFactor

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ1TQ0simple (rvanka, rvector, rrhs, domega)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8
    real(DP) :: daux
    !REAL(DP), DIMENSION(3,3) :: Dmult

    ! Local arrays for informations about one element
    real(DP), dimension(4) :: AA,BB1,BB2,DD1,DD2
    real(DP), dimension(9) :: FF,UU
    integer, dimension(4) :: idofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanka%Dmultipliers(:,:)

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    do iel=1,NEL

      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      do inode=1,4

        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))

        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANKA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let us hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the four velocity unknowns and the
        ! pressure unknown on the current element; these five unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity
        ! components and B being an (2*4) x 1 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |        |             |        |        |
        !   --|---X----|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)

          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux

          ! For support of scaled matrices, use the following lines; currently switched off.
          !FF(inode)       = FF(inode)      -Dmult(1,3)*p_DB1(ib)*daux
          !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,3)*p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)

            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,3)*p_DB1(ib)
            !BB2(inode) = Dmult(2,3)*p_DB2(ib)
            !DD1(inode) = Dmult(3,1)*p_DD1(ib)
            !DD2(inode) = Dmult(3,2)*p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,4
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ1TQ0simpleSol (rvanka, rvector, rrhs, domega, rtempVector)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! Solution-based algorithm which does not use the defect correction
  ! approach (for FEAT1-compatibility).
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(inout) :: rvector

  ! A temporary vector in the size and structure of rvector
  type(t_vectorBlock), intent(inout) :: rtempVector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8
    real(DP) :: daux
    !REAL(DP), DIMENSION(3,3) :: Dmult

    ! Local arrays for informations about one element
    real(DP), dimension(4) :: AA,BB1,BB2,DD1,DD2
    real(DP), dimension(9) :: FF,UU
    integer, dimension(4) :: idofGlobal

    ! WARNING: DOCUMENTATION PARTIALLY WRONG AND INCOMPLETE!
    ! Preconditioner was build from FEAT1 in a quick-and-dirty way...

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanka%Dmultipliers(:,:)

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Backup the solution vector
    call lsysbl_copyVector (rvector,rtempVector)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    do iel=1,NEL

      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      do inode=1,4

        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))

        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANKA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let us hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the four velocity unknowns and the
        ! pressure unknown on the current element; these five unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity
        ! components and B being an (2*4) x 1 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |        |             |        |        |
        !   --|---X----|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)

          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        end do

        ! The diagonal entry was also subtracted, but this is not
        ! desired in this approach. We therefore revert the subtraction
        ! of the diagonal element.
        FF(inode)       = FF(inode)       + AA(inode)*p_Dvector(idof)
        FF(inode+lofsv) = FF(inode+lofsv) + AA(inode)*p_Dvector(idof+ioffsetv)

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1

        do ib = ib1,ib2
          ! Subtract contributions from the RHS which do not belong to our element.
          if (p_KcolB(ib) .ne. IEL) then

            J = p_KcolB(ib)
            daux = p_Dvector(j+ioffsetp)
            FF(inode)       = FF(inode)      -p_DB1(ib)*daux
            FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux

            ! For support of scaled matrices, use the following lines; currently switched off.
            !FF(inode)       = FF(inode)      -Dmult(1,3)*p_DB1(ib)*daux
            !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,3)*p_DB2(ib)*daux

          else

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)

            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,3)*p_DB1(ib)
            !BB2(inode) = Dmult(2,3)*p_DB2(ib)
            !DD1(inode) = Dmult(3,1)*p_DD1(ib)
            !DD2(inode) = Dmult(3,2)*p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            !
            ! COMMENTED OUT, we are working solution based!
            ! FF(1+lofsp) = FF(1+lofsp) &
            !             - DD1(inode)*p_Dvector(idof) &
            !             - DD2(inode)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            ! EXIT
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,4
        p_Dvector(idofGlobal(inode)) = UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) = UU(inode+lofsv)
      end do

      p_Dvector(iel+ioffsetp) = UU(1+lofsp)

    end do ! iel

    ! The final vector is a mean between the old and the new vector.
    call lsysbl_vectorLinearComb (rtempVector,rvector,1.0_DP-domega,domega)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ1TQ0simpleConf (rvanka, rvector, rrhs, domega,IelementList)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ1TQ0simpleConf is the same as vanka_2DSPQ1TQ0simple except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8
    real(DP) :: daux

    ! Local arrays for informations about one element
    real(DP), dimension(4) :: AA,BB1,BB2,DD1,DD2
    real(DP), dimension(9) :: FF,UU
    integer, dimension(4) :: idofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      do inode=1,4

        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the four velocity unknowns and the
        ! pressure unknown on the current element; these five unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity
        ! components and B being an (2*4) x 1 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |        |             |        |        |
        !   --|---X----|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,4
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)

    end do ! iel

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine vanka_getcorr_2DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,DD1,DD2)

!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two
  ! velocity vectors and one pressure vector. It is used as auxiliary
  ! routine in the simple VANKA solver to calculate an update vector
  ! for velocity/pressure.
!</description>

!<input>
  ! Diagonal elements of the local system matrix.
  real(DP), dimension(4), intent(in) :: AA

  ! Entries in the submatrix B1.
  real(DP), dimension(4), intent(in) :: BB1

  ! Entries in the submatrix B2.
  real(DP), dimension(4), intent(in) :: BB2

  ! Entries in the submatrix D1.
  real(DP), dimension(4), intent(in) :: DD1

  ! Entries in the submatrix D2.
  real(DP), dimension(4), intent(in) :: DD2

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  real(DP), dimension(9), intent(in) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  real(DP), dimension(9), intent(out) :: UU
!</output>

!</subroutine>

    ! local variables

    integer :: inode
    real(DP) :: PP,dpres
    real(DP), dimension(9) :: AI,dff

    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA(1)                                                   BB1(1) )
    !     (        AA(2)                                            BB1(2) )
    !     (               AA(3)                                     BB1(3) )
    !     (                      AA(4)                              BB1(4) )
    !     (                             AA(1)                       BB2(1) )
    !     (                                    AA(2)                BB2(2) )
    !     (                                           AA(3)         BB2(3) )
    !     (                                                  AA(4)  BB2(4) )
    !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial!
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    do inode=1,4

      ! Quick check if everything is ok - we do not want to divide by 0.
      if (AA(inode)*AA(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      ! AI(.) saves the diagonal matrix S^-1:

      AI(inode)=1E0_DP/AA(inode)

    end do

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 4x1 matrix
    ! b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = 0.0_DP
    dff = FF

    do inode = 1,4
      dpres        = dpres &
                   - AI(inode)*(DD1(inode)*BB1(inode)+DD2(inode)*BB2(inode))
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI(inode)*(DD1(inode)*dff(inode)+DD2(inode)*dff(inode+lofsv))
    end do

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF`s are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    if (dpres*dpres .lt. 1E-20_DP)  then
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      return
    endif

    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:

    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP

    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:

    do inode=1,4
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
    end do

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine vanka_getcorr_2DSPQ1TQ0simple2 (UU,FF,AA1,AA2,BB1,BB2,DD1,DD2,di)

!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two
  ! velocity vectors and one pressure vector. It is used as auxiliary
  ! routine in the simple VANKA solver to calculate an update vector
  ! for velocity/pressure.
  !
  ! In contrast to vanka_getcorr_2DSPQ1TQ0simple, the two diagonal blocks
  ! in the matrix may be different from each other.
!</description>

!<input>
  ! Diagonal elements of the local system matrix A11
  real(DP), dimension(4), intent(in) :: AA1

  ! Diagonal elements of the local system matrix A22
  real(DP), dimension(4), intent(in) :: AA2

  ! Entries in the submatrix B1.
  real(DP), dimension(4), intent(in) :: BB1

  ! Entries in the submatrix B2.
  real(DP), dimension(4), intent(in) :: BB2

  ! Entries in the submatrix D1.
  real(DP), dimension(4), intent(in) :: DD1

  ! Entries in the submatrix D2.
  real(DP), dimension(4), intent(in) :: DD2

  ! Entry in the submatrix I (usually =0).
  ! This is the matrix in the diagonal block of the pressure, which is usually
  ! zero in saddle point problems.
  real(DP), intent(in) :: di

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  real(DP), dimension(9), intent(in) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  real(DP), dimension(9), intent(out) :: UU
!</output>

!</subroutine>

    ! local variables

    integer :: inode
    real(DP) :: PP,dpres
    real(DP), dimension(4) :: AI1,AI2
    real(DP), dimension(9) :: dff

    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA1(1)                                                  BB1(1) )
    !     (        AA1(2)                                           BB1(2) )
    !     (               AA1(3)                                    BB1(3) )
    !     (                      AA1(4)                             BB1(4) )
    !     (                             AA2(1)                      BB2(1) )
    !     (                                    AA2(2)               BB2(2) )
    !     (                                           AA2(3)        BB2(3) )
    !     (                                                  AA2(4) BB2(4) )
    !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4) di1    )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T i )
    !      ( D1  D2  I  )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t I ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !       D^t u + Ip = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !      D^t S^-1 (f-Bp) + Ip = g
    !
    !      <=>  Ip -D^t S^-1 B p  =  g - D^t S^-1 f
    !              ***********       **************
    !                 =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial!
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    do inode=1,4

      ! Quick check if everything is ok - we do not want to divide by 0.
      if (AA1(inode)*AA1(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      if (AA2(inode)*AA2(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      ! AI(.) saves the diagonal matrix S^-1:

      AI1(inode)=1E0_DP/AA1(inode)
      AI2(inode)=1E0_DP/AA2(inode)

    end do

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !      ( I - D^t S^-1 B ) p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 4x1 matrix
    ! b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           = (I - D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = di
    dff = FF

    do inode = 1,4
      dpres        = dpres &
                   - AI1(inode)*DD1(inode)*BB1(inode) &
                   - AI2(inode)*DD2(inode)*BB2(inode)
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI1(inode)*DD1(inode)*dff(inode) &
                   - AI2(inode)*DD2(inode)*dff(inode+lofsv)
    end do

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF`s are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    if (dpres*dpres .lt. 1E-20_DP)  then
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      return
    endif

    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:

    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP

    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:

    do inode=1,4
      UU(inode)       = AI1(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI2(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
    end do

  end subroutine

  ! ***************************************************************************
  ! 2D Navier-Stokes VANKA, simple diagonal version for fully coupled
  ! velocity matrix.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1 )
  !    ( A21  A22  B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ1TQ0simpleCoupConf (rvanka, rvector, rrhs, domega,IelementList)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ1TQ0simpleCoupConf is the same as vanka_2DSPQ1TQ0simpleConf,
  ! but supports fully coupled velocity submatrices.
  ! The matrices A11 and A22 must have the same structure. The matrices A12
  ! and A21 must also have the same structure. The structure of A11 and A12
  ! may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12,p_KdiagonalA12
    real(DP), dimension(:), pointer :: p_DA,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8
    real(DP) :: daux

    ! Local arrays for informations about one element
    real(DP), dimension(2*lofsv) :: AA,BB1,BB2,DD1,DD2
    real(DP), dimension(9) :: FF,UU
    integer, dimension(lofsv) :: idofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_DA22 => rvanka%p_DA22

    ! Structure of A12 is assumed to be the same as A21
    p_KcolA12 => rvanka%p_KcolA12
    p_KldA12 => rvanka%p_KldA12
    p_KdiagonalA12 => rvanka%p_KdiagonalA12
    p_DA12 => rvanka%p_DA12
    p_DA21 => rvanka%p_DA21

    ! Structure of B1 is assumed to be the same as B2
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      do inode=1,lofsv

        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+lofsv) = p_DA22(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the four velocity unknowns and the
        ! pressure unknown on the current element; these five unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity
        ! components and B being an (2*4) x 1 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |        |             |        |        |
        !   --|---X----|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)
        end do

        ! Tackle 'offdiagonal' matrices A12 and A21

        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_2DSPQ1TQ0CPsimple (UU,FF,AA,BB1,BB2,DD1,DD2)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,lofsv
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)

    end do ! iel

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine vanka_getcorr_2DSPQ1TQ0CPsimple (UU,FF,AA,BB1,BB2,DD1,DD2)

!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two
  ! velocity vectors and one pressure vector. It is used as auxiliary
  ! routine in the simple VANKA solver to calculate an update vector
  ! for velocity/pressure for system where the velocity is fully coupled.
!</description>

!<input>
  ! Diagonal elements of the local system matrix. (Blocks A11 and A22)
  real(DP), dimension(2*4), intent(in) :: AA

  ! Entries in the submatrix B1.
  real(DP), dimension(4), intent(in) :: BB1

  ! Entries in the submatrix B2.
  real(DP), dimension(4), intent(in) :: BB2

  ! Entries in the submatrix D1.
  real(DP), dimension(4), intent(in) :: DD1

  ! Entries in the submatrix D2.
  real(DP), dimension(4), intent(in) :: DD2

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  real(DP), dimension(9), intent(in) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  real(DP), dimension(9), intent(out) :: UU
!</output>

!</subroutine>

    ! local variables

    integer :: inode
    real(DP) :: PP,dpres
    real(DP), dimension(9) :: AI,dff

    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA(1)                                                   BB1(1) )
    !     (        AA(2)                                            BB1(2) )
    !     (               AA(3)                                     BB1(3) )
    !     (                      AA(4)                              BB1(4) )
    !     (                             AA(5)                       BB2(1) )
    !     (                                    AA(6)                BB2(2) )
    !     (                                           AA(7)         BB2(3) )
    !     (                                                  AA(8)  BB2(4) )
    !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial!
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    do inode=1,2*lofsv

      ! Quick check if everything is ok - we do not want to divide by 0.
      if (AA(inode)*AA(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      ! AI(.) saves the diagonal matrix S^-1:

      AI(inode)=1E0_DP/AA(inode)

    end do

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 4x1 matrix
    ! b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = 0.0_DP
    dff = FF

    do inode = 1,lofsv
      dpres        = dpres &
                   - AI(inode)  *(DD1(inode)*BB1(inode)) &
                   - AI(inode+4)*(DD2(inode)*BB2(inode))
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI(inode)  *(DD1(inode)*dff(inode)) &
                   - AI(inode+4)*(DD2(inode)*dff(inode+lofsv))
    end do

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF`s are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    if (dpres*dpres .lt. 1E-20_DP)  then
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      return
    endif

    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:

    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP

    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:

    do inode=1,lofsv
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode+4)*(dff(inode+lofsv)-BB2(inode)*PP)
    end do

  end subroutine

  ! ***************************************************************************
  ! 2D Navier-Stokes VANKA, 'full' version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ1TQ0fullConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ2QP1fullConf is the same as vanka_2DSPQ2QP1full except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 4      ! Q1T = 4 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! QQ0 = 1 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          4    P    2
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = edge.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 8 x 2*(NMT) matrix for the two velocity
        ! components and B~ being an (2*4) x 2 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! 2*1 pressure elements on the adjacent cells, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 2 submatrices (originally, every velocity couples with
        ! the pressure DOF`s on that cell, so we have
        ! 1 column in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            AA(inode,      2*nnvel+1) = p_DB1(ib)
            AA(inode+nnvel,2*nnvel+1) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(2*nnvel+1,inode)       = p_DD1(ib)
            AA(2*nnvel+1,inode+nnvel) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(2*nnvel+1,inode)*p_Dvector(idof) &
                        - AA(2*nnvel+1,inode+nnvel)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  .............. :: :: :: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(4,4) :: :: :: )
      !     (                            AA( 5, 5) .............. :: :: :: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA( 8, 8) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPQ1TQ0fullConf')

      end if
      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************
  ! 2D Navier-Stokes VANKA, 'full' version for fully coupled systems.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1 )
  !    ( A21  A22  B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ1TQ0fullCoupConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ1TQ0fullCoupConf is the same as vanka_2DSPQ1TQ0fullConf,
  ! but supports fully coupled velocity submatrices.
  ! The matrices A11 and A22 must have the same structure. The matrices A12
  ! and A21 must also have the same structure. The structure of A11 and A12
  ! may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12,p_KdiagonalA12
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 4      ! Q1T = 4 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! QQ0 = 1 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_DA22 => rvanka%p_DA22

    ! Structure of A12 is assumed to be the same as A21
    p_KcolA12 => rvanka%p_KcolA12
    p_KldA12 => rvanka%p_KldA12
    p_KdiagonalA12 => rvanka%p_KdiagonalA12
    p_DA12 => rvanka%p_DA12
    p_DA21 => rvanka%p_DA21

    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          4    P    2
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = edge.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:nnvel) = p_IedgesAtElement(1:nnvel,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 8 x 2*(NMT) matrix for the two velocity
        ! components and B~ being an (2*4) x 2 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! 2*1 pressure elements on the adjacent cells, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with four square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 2 submatrices (originally, every velocity couples with
        ! the pressure DOF`s on that cell, so we have
        ! 1 column in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = p_DA(ia)
              AA (inode+nnvel,k+nnvel) = p_DA22(ia)
              exit
            end if
          end do
        end do

        ! Process the 'off-diagonal' matrices A12 and A21

        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k+nnvel) = p_DA12(ia)
              AA (inode+nnvel,k) = p_DA21(ia)
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            AA(inode,      2*nnvel+1) = p_DB1(ib)
            AA(inode+nnvel,2*nnvel+1) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(2*nnvel+1,inode)       = p_DD1(ib)
            AA(2*nnvel+1,inode+nnvel) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(2*nnvel+1,inode)*p_Dvector(idof) &
                        - AA(2*nnvel+1,inode+nnvel)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  ..............    AA( 1, 5) .............. :: :: :: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B1): )
      !     (    ................ AA(4,4)    ............... AA( 4, 8) :: :: :: )
      !     ( AA(5,1)  ..............    AA( 5, 5) .............. :: :: :: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B2): )
      !     (    ................ AA(8,4)    ............... AA( 8, 8) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)

      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPQ1TQ0fullCoupConf')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************
  ! 2D Navier-Stokes VANKA, simple diagonal and full version.
  ! Supports only Q2/QP1 discretisation.
  ! Matrix must be of the form
  !
  !    ( A         B1 )
  !    (      A    B2 )
  !    ( D1^T D2^T 0  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ2QP1simple (rvanka, rvector, rrhs, domega)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer :: NVT
    integer :: NMT
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 9      ! Q2 = 9 DOF`s per velocity
    integer, parameter :: nnpressure = 3 ! QP1 = 3 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF
    real(DP) :: dmult11, dmult22, dmult13, dmult31

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j,isubdof
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    dmult11 = rvanka%Dmultipliers(1,1)
    dmult22 = rvanka%Dmultipliers(2,2)
    dmult31 = rvanka%Dmultipliers(3,1)
    dmult13 = rvanka%Dmultipliers(1,3)

    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    do iel=1,NEL

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          8    9    6
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Put on AA(.) the diagonal entry of matrix A -- the 1st and the
        ! 2nd block
        AA(inode,inode) = dmult11*p_DA(p_KdiagonalA(idof))
        AA(inode+nnvel,inode+nnvel) = dmult22*p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity
        ! components and B~ being an (2*9) x 12 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most
        ! 4*3 pressure elements on the neighbour cell, so we have
        ! 12 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! or
        !
        !              IEL
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF`s on that cell, so we have
        ! 3 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp) * dmult13
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of the original B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF`s), 6 (if the velocity DOF is an edge with
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on
        ! the boundary and there is no neighbour, or if it is the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.

        do ib = ib1,ib2

          if (p_KcolB(ib) .eq. IEL) then
            isubdof = 1
          else if (p_KcolB(ib) .eq. IEL+NEL) then
            isubdof = 2
          else if (p_KcolB(ib) .eq. IEL+NEL*2) then
            isubdof = 3
          else
            ! Cycle the loop - the entry belongs to another
            ! element, not to the current one
            cycle
          end if

          J = p_KcolB(ib)

          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = dmult13*p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = dmult13*p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = dmult31*p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = dmult31*p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1) :: :: :: )
      !     (          ..                                               :AA :: )
      !     (               ..                                          :(B1): )
      !     (                     AA(9,9) :: :: :: )
      !     (                            AA(10,10) :: :: :: )
      !     (                                      ..                   :AA :: )
      !     (                                           ..              :(B2): )
      !     (                                                AA(18,18) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPQ2QP1simple')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ2QP1simpleConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ2QP1simpleConf is the same as vanka_2DSPQ2QP1simple except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof
    logical :: bsuccess

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer :: NVT
    integer :: NMT
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 9      ! Q2 = 9 DOF`s per velocity
    integer, parameter :: nnpressure = 3 ! QP1 = 3 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(2*nnvel) :: AA
    real(DP), dimension(nnpressure,2*nnvel) :: DD
    real(DP), dimension(2*nnvel,nnpressure) :: BB
    real(DP), dimension(nnpressure) :: CC
    real(DP), dimension(nnld) :: FF,UU

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j,isubdof
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    ! Currently no C supported.
    CC(:) = 0

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the B and D-matrices.
      ! Note: THis is not necessary as all entries of BB and DD
      ! will be overwritten later! Reason: Because of the FEM
      ! coupling, all entries are affected.
      ! BB(:,:) = 0.0_DP
      ! DD(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          8    9    6
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Put on AA(.) the diagonal entry of matrix A -- the 1st and the
        ! 2nd block
        AA(inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+nnvel) = p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity
        ! components and B~ being an (2*9) x 12 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most
        ! 4*3 pressure elements on the neighbour cell, so we have
        ! 12 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! or
        !
        !              IEL
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF`s on that cell, so we have
        ! 3 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA(ia)*p_Dvector(J+ioffsetv)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of the original B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF`s), 6 (if the velocity DOF is an edge with
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on
        ! the boundary and there is no neighbour, or if it is the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.

        do ib = ib1,ib2

          if (p_KcolB(ib) .eq. IEL) then
            isubdof = 1
          else if (p_KcolB(ib) .eq. IEL+NEL) then
            isubdof = 2
          else if (p_KcolB(ib) .eq. IEL+NEL*2) then
            isubdof = 3
          else
            ! Cycle the loop - the entry belongs to another
            ! element, not to the current one
            cycle
          end if

          J = p_KcolB(ib)

          ! Get the entries in the B-matrices
          BB(inode,      isubdof) = p_DB1(ib)
          BB(inode+nnvel,isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          DD(isubdof,inode)       = p_DD1(ib)
          DD(isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - DD(isubdof,inode)*p_Dvector(idof) &
                      - DD(isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1) :: :: :: )
      !     (          ..                                               :AA :: )
      !     (               ..                                          :(B1): )
      !     (                     AA(9,9) :: :: :: )
      !     (                            AA(10,10) :: :: :: )
      !     (                                      ..                   :AA :: )
      !     (                                           ..              :(B2): )
      !     (                                                AA(18,18) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      !call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      call qsol_solveDiagSchurComp (2*nnvel,nnpressure,UU,FF,AA,BB,DD,CC,bsuccess)

      if (bsuccess) then

        ! Ok, we got the update vector in UU. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * UU(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * UU(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * UU(3+lofsp)

      end if

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ2QP1simpleCoupConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ2QP1simpleCoupConf is the same as vanka_2DSPQ2QP1simpleConf but
  ! allows different matrices in A12/A21.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof
    logical :: bsuccess

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA,p_DA22
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12,p_KdiagonalA12
    real(DP), dimension(:), pointer :: p_DA12,p_DA21
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer :: NVT
    integer :: NMT
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 9      ! Q2 = 9 DOF`s per velocity
    integer, parameter :: nnpressure = 3 ! QP1 = 3 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(2*nnvel) :: AA
    real(DP), dimension(nnpressure,2*nnvel) :: DD
    real(DP), dimension(2*nnvel,nnpressure) :: BB
    real(DP), dimension(nnpressure) :: CC
    real(DP), dimension(nnld) :: FF,UU

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j,isubdof
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_DA22 => rvanka%p_DA22
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    p_KcolA12 => rvanka%p_KcolA12
    p_KldA12 => rvanka%p_KldA12
    p_KdiagonalA12 => rvanka%p_KdiagonalA12
    p_DA12 => rvanka%p_DA12
    p_DA21 => rvanka%p_DA21

    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    ! Currently no C supported.
    CC(:) = 0

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the B and D-matrices.
      ! Note: THis is not necessary as all entries of BB and DD
      ! will be overwritten later! Reason: Because of the FEM
      ! coupling, all entries are affected.
      ! BB(:,:) = 0.0_DP
      ! DD(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          8    9    6
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Put on AA(.) the diagonal entry of matrix A -- the 1st and the
        ! 2nd block
        AA(inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+nnvel) = p_DA22(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity
        ! components and B~ being an (2*9) x 12 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most
        ! 4*3 pressure elements on the neighbour cell, so we have
        ! 12 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! or
        !
        !              IEL
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF`s on that cell, so we have
        ! 3 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)
        end do

        ! Tackle 'offdiagonal' matrices A12 and A21

        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of the original B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF`s), 6 (if the velocity DOF is an edge with
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on
        ! the boundary and there is no neighbour, or if it is the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.

        do ib = ib1,ib2

          if (p_KcolB(ib) .eq. IEL) then
            isubdof = 1
          else if (p_KcolB(ib) .eq. IEL+NEL) then
            isubdof = 2
          else if (p_KcolB(ib) .eq. IEL+NEL*2) then
            isubdof = 3
          else
            ! Cycle the loop - the entry belongs to another
            ! element, not to the current one
            cycle
          end if

          J = p_KcolB(ib)

          ! Get the entries in the B-matrices
          BB(inode,      isubdof) = p_DB1(ib)
          BB(inode+nnvel,isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          DD(isubdof,inode)       = p_DD1(ib)
          DD(isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - DD(isubdof,inode)*p_Dvector(idof) &
                      - DD(isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1) :: :: :: )
      !     (          ..                                               :AA :: )
      !     (               ..                                          :(B1): )
      !     (                     AA(9,9) :: :: :: )
      !     (                            AA(10,10) :: :: :: )
      !     (                                      ..                   :AA :: )
      !     (                                           ..              :(B2): )
      !     (                                                AA(18,18) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      !call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
      call qsol_solveDiagSchurComp (2*nnvel,nnpressure,UU,FF,AA,BB,DD,CC,bsuccess)

      if (bsuccess) then

        ! Ok, we got the update vector in UU. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * UU(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * UU(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * UU(3+lofsp)
      end if

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ2QP1full (rvanka, rvector, rrhs, domega)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer :: NVT
    integer :: NMT
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 9      ! Q2 = 9 DOF`s per velocity
    integer, parameter :: nnpressure = 3 ! QP1 = 3 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,isubdof,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements

    do iel=1,NEL

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          8    9    6
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity
        ! components and B~ being an (2*9) x 12 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most
        ! 4*3 pressure elements on the adjacent cells, so we have
        ! 12 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! or
        !
        !              IEL
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF`s on that cell, so we have
        ! 3 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF`s), 6 (if the velocity DOF is an edge with
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on
        ! the boundary and there is no neighbour, or if it is the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.

        do ib = ib1,ib2

          if (p_KcolB(ib) .eq. IEL) then
            isubdof = 1
          else if (p_KcolB(ib) .eq. IEL+NEL) then
            isubdof = 2
          else if (p_KcolB(ib) .eq. IEL+NEL*2) then
            isubdof = 3
          else
            ! Cycle the loop - the entry belongs to another
            ! element, not to the current one
            cycle
          end if

          J = p_KcolB(ib)

          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  .............. :: :: :: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(9,9) :: :: :: )
      !     (                            AA(10,10) .............. :: :: :: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA(18,18) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPQ2QP1full')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ2QP1fullConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ2QP1fullConf is the same as vanka_2DSPQ2QP1full except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer :: NVT
    integer :: NMT
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 9      ! Q2 = 9 DOF`s per velocity
    integer, parameter :: nnpressure = 3 ! QP1 = 3 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,isubdof,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          8    9    6
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity
        ! components and B~ being an (2*9) x 12 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most
        ! 4*3 pressure elements on the adjacent cells, so we have
        ! 12 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! or
        !
        !              IEL
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF`s on that cell, so we have
        ! 3 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF`s), 6 (if the velocity DOF is an edge with
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on
        ! the boundary and there is no neighbour, or if it is the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.

        do ib = ib1,ib2

          if (p_KcolB(ib) .eq. IEL) then
            isubdof = 1
          else if (p_KcolB(ib) .eq. IEL+NEL) then
            isubdof = 2
          else if (p_KcolB(ib) .eq. IEL+NEL*2) then
            isubdof = 3
          else
            ! Cycle the loop - the entry belongs to another
            ! element, not to the current one
            cycle
          end if

          J = p_KcolB(ib)

          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  .............. :: :: :: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(9,9) :: :: :: )
      !     (                            AA(10,10) .............. :: :: :: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA(18,18) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPQ2QP1fullConf')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPQ2QP1fullCoupConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier-Stokes problems with Q2/Q1 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ2QP1fullCoupConf is the same as vanka_2DSPQ2QP1fullCoupConf
  ! but also allows for different matrices A11,A12,A21 and A22.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA11,p_KcolA12
    integer, dimension(:), pointer :: p_KldA11
    integer, dimension(:), pointer :: p_KldA12
    real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer :: NVT
    integer :: NMT
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 9      ! Q2 = 9 DOF`s per velocity
    integer, parameter :: nnpressure = 3 ! QP1 = 3 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,isubdof,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(dp) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA11 => rvanka%p_KcolA
    p_KcolA12 => rvanka%p_KcolA12
    p_KldA11 => rvanka%p_KldA
    p_KldA12 => rvanka%p_KldA12
    p_DA11 => rvanka%p_DA
    p_DA12 => rvanka%p_DA12
    p_DA21 => rvanka%p_DA21
    p_DA22 => rvanka%p_DA22
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NVT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NVT
    NMT = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NMT
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !      old  old  old            new  new  new
    !        X---X---X                X---X---X
    !        |       |                |       |
    !    old X   X   X old   -->  new X   X   X new
    !        |   1   |                |   1   |
    !        X---X---X                X---X---X
    !      old  old  old            new  new  new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !      new  new new old  old         new  new newer new new
    !        X---X---X---X---X             X---X---|---X---X
    !        |       |       |             |   1   |   2   |
    !    new X   X   X   X   X old --> new X   X   X   X   X new
    !        |   1   |new 1  |             |       |newer  |
    !        X---X---X---X---X             X---X---X---X---X
    !      new  new new old  old         new  new newer new new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       4----7----3
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          8    9    6
      ! |         |                       |    P1-3 |
      ! +---------+                       1----5----2
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
      FF(2+lofsp) = p_Drhs(iel+NEL+ioffsetp)
      FF(3+lofsp) = p_Drhs(iel+2*NEL+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = corner vertex, DOF 5..8 = edge, DOF 9 = element.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IverticesAtElement(1:4,iel)
      IdofGlobal(5:8) = p_IedgesAtElement(1:4,iel)+NVT
      IdofGlobal(9)   = NVT+NMT+iel

      ! Loop over all 9 U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 18 x 2*(NVT+NMT+NEL) matrix for the two velocity
        ! components and B~ being an (2*9) x 12 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 18 x 12 matrix: As every velocity couples with at most
        ! 4*3 pressure elements on the adjacent cells, so we have
        ! 12 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! or
        !
        !              IEL
        ! |--------|--------|
        ! |        |        |
        ! |   Q1   |   P    |
        ! |        |        |
        ! |--------X--------|   or X a vertex or an edge on the boundary.
        ! |        |        |
        ! |   Q2   |   Q3   |
        ! |        |        |
        ! |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 9 x 9. The 18 x 12-matrix B~ reduces to
        ! two 9 x 3 submatrices (originally, every velocity couples with
        ! the 3 pressure DOF`s on that cell, so we have
        ! 3 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui
        !
        ! Diagonal matrices
        ia1 = p_KldA11(idof)
        ia2 = p_KldA11(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA11(ia)
          FF(inode)       = FF(inode)      -p_DA11(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = p_DA11(ia)
              AA (inode+nnvel,k+nnvel) = p_DA22(ia)
              exit
            end if
          end do
        end do

        ! Offdiagonal matrices
        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k+nnvel) = p_DA12(ia)
              AA (inode+nnvel,k) = p_DA21(ia)
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! In the next loop we have to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most 12 entries:
        !
        !      IEL                              IEL
        !   |-------|              +--------X-------+
        !   |       |              |        |       |
        !   |   P1  |       or     |   P2   X   X   |
        !   |       |              |        |   P1  |
        ! --X---X---X--            +--------X-------+
        !                          |        |       |
        !                          |   P3   |   P4  |
        !                          |        |       |
        !                          +--------+-------+
        !
        ! Either 12 (for corner DOF`s), 6 (if the velocity DOF is an edge with
        ! two neighbouring elements) or 3 (if the velocity DOF is at an edge on
        ! the boundary and there is no neighbour, or if it is the element midpoint).
        !
        ! 3 of these 12 entries in each line come into our 'local' B-matrices.

        do ib = ib1,ib2

          if (p_KcolB(ib) .eq. IEL) then
            isubdof = 1
          else if (p_KcolB(ib) .eq. IEL+NEL) then
            isubdof = 2
          else if (p_KcolB(ib) .eq. IEL+NEL*2) then
            isubdof = 3
          else
            ! Cycle the loop - the entry belongs to another
            ! element, not to the current one
            cycle
          end if

          J = p_KcolB(ib)

          ! Get the entries in the B-matrices
          AA(inode,      2*nnvel+isubdof) = p_DB1(ib)
          AA(inode+nnvel,2*nnvel+isubdof) = p_DB2(ib)

          ! The same way, get DD1 and DD2.
          ! Note that DDi has exacty the same matrix structrure as BBi and is noted
          ! as 'transposed matrix' only because of the transposed-flag.
          ! So we can use "ib" as index here to access the entry of DDi:
          AA(2*nnvel+isubdof,inode)       = p_DD1(ib)
          AA(2*nnvel+isubdof,inode+nnvel) = p_DD2(ib)

          ! Build the pressure entry in the local defect vector:
          !   f_i = (f_i-Aui) - D_i pi
          ! or more precisely (as D is roughly B^T):
          !   f_i = (f_i-Aui) - (B^T)_i pi
          FF(isubdof+lofsp) = FF(isubdof+lofsp) &
                      - AA(2*nnvel+isubdof,inode)*p_Dvector(idof) &
                      - AA(2*nnvel+isubdof,inode+nnvel)*p_Dvector(idof+ioffsetv)
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  ..............    AA(1,10) ............... :: :: :: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B1): )
      !     (    ................ AA(9,9)    ............... AA(10,18) :: :: :: )
      !     ( AA(10,1) ..............    AA(10,10) .............. :: :: :: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B2): )
      !     (    ............... AA(18,9)    ............... AA(18,18) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
        p_Dvector(iel+NEL+ioffsetp) = p_Dvector(iel+NEL+ioffsetp) + &
                                      domega * FF(2+lofsp)
        p_Dvector(iel+2*NEL+ioffsetp) = p_Dvector(iel+2*NEL+ioffsetp) + &
                                        domega * FF(3+lofsp)
      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPQ2QP1fullConf')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************
  ! 2D VANKA, 'full' version for fully coupled Navier-Stokes.
  ! Extended version with support for diagonal matrices in the pressure and
  ! arbitrary scaled submatrices.
  ! Supports only Q1~/Q0.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1  )
  !    ( A21  A22  B2  )
  !    ( D1^T D2^T I1  )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  !
  ! I1 is a diagonal matrix in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! do not exist, what results in a '0' block in these positions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DNSQ1TQ0fullCoupConfExt (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier Stokes optimal control problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DNSSQ1TQ0fullCoupConf supports fully coupled velocity submatrices.
  ! The matrices A11, A22 must have the same structure.
  ! The matrices A12 and A21 must have the same structure.
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA11
    integer, dimension(:), pointer :: p_KldA11
    real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2
    real(DP), dimension(:), pointer :: p_Da33
    integer, dimension(:), pointer :: p_KdiagonalA33

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 4      ! Q1T = 4 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! QQ0 = 1 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel + nnpressure
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! Offsets of the 'local' solution parts in the 'local' solution vector
    integer, parameter :: lofsu = 0
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! Offset information in arrays.
    integer :: ioffsetu,ioffsetv,ioffsetp,j

    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA11 => rvanka%p_KcolA
    p_KldA11 => rvanka%p_KldA
    p_DA11 => rvanka%p_DA
    p_DA22 => rvanka%p_DA22
    if (.not. associated(p_DA22)) p_DA22 => p_DA11

    ! Structure of A12 is assumed to be the same as A21.
    ! Get A12 and A21 -- except for if the multipliers are =0, then
    ! we switch them off by nullifying the pointers.
    if (rvanka%Dmultipliers(1,2) .ne. 0.0_DP) then
      p_KcolA12 => rvanka%p_KcolA12
      p_KldA12 => rvanka%p_KldA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    else
      nullify(p_KcolA12)
      nullify(p_KldA12)
      nullify(p_DA12 )
      nullify(p_DA21 )
    end if

    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Diagonal submatrices A33 and A66 (if they exist)
    if (rvanka%Dmultipliers(3,3) .ne. 0.0_DP) then
      p_Da33 => rvanka%p_DA33
      p_KdiagonalA33 => rvanka%p_KdiagonalA33
    else
      nullify(p_Da33)
      nullify(p_KdiagonalA33)
    end if

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          4    P    2
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !
      !
      ! Fetch the pressure P on the current element into FF.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      ! Get the pressure
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = edge.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode.

        ! Primal equation
        FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u)  = (f  )
        !    [ B^t 0 ] (p)    (g  )
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in rectangular
        ! systems of the form
        !
        !    [ === A^ === B~  ] (| ) = (f1 )
        !    [ B~^t       I1~ ] (u )   (f2 )
        !                       (| )   (g  )
        !                       (p )
        !
        !
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! 2*1 pressure elements on the adjacent cells, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL.
        !
        !  (d1 ) = (f1 ) -  [ === A^ === B~  ] (| )
        !  (d2 )   (f2 )    [ B~^t       I1~ ] (u )
        !  (dg )   (g  )                       (| )
        !                                      (p )
        !
        ! Extract those entries in the A-, B- and M-matrices to our local
        ! matrix AA, which belong to DOF`s in our current solution vector.
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA11(idof)
        ia2 = p_KldA11(idof+1)-1
        do ia = ia1,ia2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - ( A11  .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .  A22  .   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   ) ( p  )

          J = p_KcolA11(ia)

          ! Primal equation:
          FF(inode+lofsu) = FF(inode+lofsu) &
                          - rvanka%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
          FF(inode+lofsv) = FF(inode+lofsv) &
                          - rvanka%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanka%Dmultipliers(1,1)
              AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanka%Dmultipliers(2,2)
              exit
            end if
          end do
        end do

        ! Process the 'off-diagonal' matrices A12 and A21

        if (associated(p_KldA12)) then
          ia1 = p_KldA12(idof)
          ia2 = p_KldA12(idof+1)-1
          do ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .  A12  .  ) ( u  )
            !   ( dv  )   ( dv  )   ( A21  .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .  ) ( p  )

            J = p_KcolA12(ia)
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanka%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanka%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanka%Dmultipliers(1,2)
                AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanka%Dmultipliers(2,1)
                exit
              end if
            end do
          end do
        end if

        ! Process A33 if it exists

        if (associated(p_KdiagonalA33)) then

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   I1  ) ( p  )
          !
          ! IEL is the pressure DOF which we have to tackle.

          daux = rvanka%Dmultipliers(3,3)
          FF(1+lofsp) = FF(1+lofsp) &
                      - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
          AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

        end if

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .  B1   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .  B2   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   ) ( p  )

          J = p_KcolB(ib)

          daux = p_Dvector(j+ioffsetp)
          FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanka%Dmultipliers(1,3)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanka%Dmultipliers(2,3)

          ! Do not incorporate the B-matrices into AA yet; this will come later!
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   ) ( v  )
          !   ( dp  )   ( dp  )   ( D1  D2   .   ) ( p  )
          !
          ! In AA, we simultaneously set up (locally):
          !
          !   (  .   .  B1   )
          !   (  .   .  B2   )
          !   ( D1  D2   .   )

          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices.
            ! Primal equation
            AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanka%Dmultipliers(1,3)
            AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanka%Dmultipliers(2,3)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanka%Dmultipliers(3,1)
            AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanka%Dmultipliers(3,2)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                        - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! For C, we use our local AA, i.e. applying C^{-1} means to
      ! solve the local system AA dd = FF for dd. The local defect dd is then
      ! added back to the global solution vector.

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y

        do inode=1,nnvel
          ! Update of the primal velocity vectors
          p_Dvector(idofGlobal(inode)+ioffsetu) &
            = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)

      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNSSQ1TQ0fullCoupConf')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

! *****************************************************************************
! Problem class: VANKA for STEADY OPTIMAL CONTROL PROBLEMS
! *****************************************************************************

!<subroutine>

  subroutine vanka_init2DNavierStokesOptC (rmatrix,rvanka)

!<description>
  ! Initialises the VANKA variant for 2D Navier-Stokes problems
  ! for conformal discretisations.
  ! Checks if the "2D-Navier-Stokes" VANKA variant
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! The substructure rvanka%rvanka2DNavSt is intitialised according
  ! to the information provided in rmatrix.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in), target :: rmatrix
  !</input>

!<inputoutput>
  ! t_vankaPointer2DSPNavSt structure that saves algorithm-specific parameters.
  type(t_vanka), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

    integer :: i,j
    type(t_blockDiscretisation), pointer :: p_rblockDiscr

    ! Matrix must be 6x6.
    if ((rmatrix%nblocksPerCol .ne. 6) .or. (rmatrix%nblocksPerRow .ne. 6)) then
      call output_line ('System matrix is not 6x6.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    ! A(1:2,1:3), A(4:5,4:6) must not be virtually transposed and of format 9.
    ! A(3,:),A(6,:) must be (virtually) transposed. All matrices must be double precision.
    do i=1,6
      do j=1,6

        if (lsysbl_isSubmatrixPresent (rmatrix,i,j)) then

          if ( ((i .ge. 1) .and. (i .le. 2)) .or. &
               ((i .ge. 4) .and. (i .le. 5)) ) then
            if (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .ne. 0) then
              call output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
              call sys_halt()
            end if
          else
            if ((i .ne. j) .and. &
                (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                 .eq. 0)) then
              call output_line ('B1/B2 submatrices must be virtually',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
              call output_line ('transposed (LSYSSC_MSPEC_TRANSPOSED)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
              call sys_halt()
            end if
          end if

          if ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX9)) then
            call output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
            call sys_halt()
          end if

          if (rmatrix%RmatrixBlock(i,j)%cdataType .ne. ST_DOUBLE) then
            call output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
            call sys_halt()
          end if

        end if ! neq != 0
      end do
    end do

    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,3)%NA .ne. rmatrix%RmatrixBlock(3,1)%NA) .or. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .ne. rmatrix%RmatrixBlock(3,1)%NCOLS)) then
      call output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,3)%NA .ne. rmatrix%RmatrixBlock(3,2)%NA) .or. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .ne. rmatrix%RmatrixBlock(3,2)%NCOLS)) then
      call output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    ! The structure of A(4,6) must be identical to A(6,4) and
    ! that of A(5,6) must be identical to A(6,5).
    if ((rmatrix%RmatrixBlock(4,6)%NA .ne. rmatrix%RmatrixBlock(6,4)%NA) .or. &
        (rmatrix%RmatrixBlock(4,6)%NEQ .ne. rmatrix%RmatrixBlock(6,4)%NCOLS)) then
      call output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    ! Primal and dual B- and D-matrices must share the same data.
    ! The matrices may be switched off by dscaleFactor=0, but they must have the
    ! same data arrays! Just for 'efficiency'...
    if ((rmatrix%RmatrixBlock(5,6)%NA .ne. rmatrix%RmatrixBlock(6,5)%NA) .or. &
        (rmatrix%RmatrixBlock(5,6)%NEQ .ne. rmatrix%RmatrixBlock(6,5)%NCOLS)) then
      call output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    if (.not. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(1,3),rmatrix%RmatrixBlock(4,6))) then
      call output_line ('Content of primal and dual B1-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    if (.not. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(2,3),rmatrix%RmatrixBlock(5,6))) then
      call output_line ('Content of primal and dual B2-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    if (.not. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(3,1),rmatrix%RmatrixBlock(6,4))) then
      call output_line ('Content of primal and dual D1-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    if (.not. lsyssc_isMatrixContentShared( &
        rmatrix%RmatrixBlock(3,2),rmatrix%RmatrixBlock(6,5))) then
      call output_line ('Content of primal and dual D2-matrix not shared!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    ! Fill the output structure with data of the matrices.
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),&
        rvanka%rvanka2DNavStOptC%p_DB1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),&
        rvanka%rvanka2DNavStOptC%p_DB2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),&
        rvanka%rvanka2DNavStOptC%p_DD1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),&
        rvanka%rvanka2DNavStOptC%p_DD2)
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),&
        rvanka%rvanka2DNavStOptC%p_KcolB)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), &
        rvanka%rvanka2DNavStOptC%p_KldB )
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanka%rvanka2DNavStOptC%p_KcolA11)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanka%rvanka2DNavStOptC%p_KldA11 )
    if (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanka%rvanka2DNavStOptC%p_KdiagonalA11)
    else
      rvanka%rvanka2DNavStOptC%p_KdiagonalA11 => rvanka%rvanka2DNavStOptC%p_KldA11
    end if

    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanka%rvanka2DNavStOptC%p_DA11 )

    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
        rvanka%rvanka2DNavStOptC%p_DA22 )

    ! What is with A12 and A21? Do they exist? With a scale factor = 1.0?
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
          rvanka%rvanka2DNavStOptC%p_DA12 )

      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,2),&
          rvanka%rvanka2DNavStOptC%p_KcolA12)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,2), &
          rvanka%rvanka2DNavStOptC%p_KldA12 )

      ! Get the structure. It is assumed that A12 and A21 have the same!
      if (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,2), &
                                rvanka%rvanka2DNavStOptC%p_KdiagonalA12)
      else
        rvanka%rvanka2DNavStOptC%p_KdiagonalA12 => rvanka%rvanka2DNavStOptC%p_KldA12
      end if

      if (.not. lsysbl_isSubmatrixPresent(rmatrix,2,1)) then
        call output_line ('If A12 is given, A21 must also be given!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
        call sys_halt()
      end if

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
          rvanka%rvanka2DNavStOptC%p_DA21 )
    end if

    ! Now we come to the dual equation -> A(4:6,4:6).
    ! We assume that A44 and A55 has the same structure as A11 and A22.

    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,4),&
        rvanka%rvanka2DNavStOptC%p_DA44 )

    ! What is with A55? Is it the same as A44?
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(5,5),&
        rvanka%rvanka2DNavStOptC%p_DA55 )

    ! What is with A12 and A21? Do they exist? With a scale factor != 0.0?
    if (lsysbl_isSubmatrixPresent(rmatrix,4,5)) then

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,5),&
          rvanka%rvanka2DNavStOptC%p_DA45 )

      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(4,5),&
          rvanka%rvanka2DNavStOptC%p_KcolA45)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,5), &
          rvanka%rvanka2DNavStOptC%p_KldA45 )

      ! Get the structure. It is assumed that A12 and A21 have the same!
      if (rmatrix%RmatrixBlock(4,5)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(4,5), &
                                rvanka%rvanka2DNavStOptC%p_KdiagonalA45)
      else
        rvanka%rvanka2DNavStOptC%p_KdiagonalA45 => rvanka%rvanka2DNavStOptC%p_KldA45
      end if

      if (.not. lsysbl_isSubmatrixPresent(rmatrix,5,4)) then
        call output_line ('If A45 is given, A54 must also be given!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
        call sys_halt()
      end if

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(5,4),&
          rvanka%rvanka2DNavStOptC%p_DA54 )
    end if

    ! Is there an identity matrix present at A(3,3) and/or A(6,6)?
    if (lsysbl_isSubmatrixPresent(rmatrix,3,3)) then
      ! The matrix must be of format 9.
      if (rmatrix%RmatrixBlock(3,3)%dscaleFactor .ne. 0.0_DP) then
        call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
            rvanka%rvanka2DNavStOptC%p_DA33 )

        if (rmatrix%RmatrixBlock(3,3)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
          call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(3,3), &
                                  rvanka%rvanka2DNavStOptC%p_KdiagonalA33)
        else
          call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,3), &
                                  rvanka%rvanka2DNavStOptC%p_KdiagonalA33)
        end if

      end if
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,6,6)) then
      if (rmatrix%RmatrixBlock(6,6)%dscaleFactor .ne. 0.0_DP) then
        call lsyssc_getbase_double(rmatrix%RmatrixBlock(6,6),&
            rvanka%rvanka2DNavStOptC%p_DA66 )

        if (rmatrix%RmatrixBlock(6,6)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
          call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(6,6), &
                                  rvanka%rvanka2DNavStOptC%p_KdiagonalA66)
        else
          call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(6,6), &
                                  rvanka%rvanka2DNavStOptC%p_KdiagonalA66)
        end if
      end if
    end if

    ! Get the mass matrix/matrices -- if they are present.
    ! It is assumed that all mass matrices are the same except for their
    ! multiplication factors!
    if (lsysbl_isSubmatrixPresent (rmatrix,1,4)) then
      if (rmatrix%RmatrixBlock(1,4)%cmatrixFormat .eq. LSYSSC_MATRIXD) then
        call output_line ('Lumped mass matrices not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
        call sys_halt()
      end if

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,4),&
          rvanka%rvanka2DNavStOptC%p_DM14 )

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,5),&
          rvanka%rvanka2DNavStOptC%p_DM25 )

      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,4),&
          rvanka%rvanka2DNavStOptC%p_KcolM)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,4), &
          rvanka%rvanka2DNavStOptC%p_KldM )

      if (rmatrix%RmatrixBlock(1,4)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,4), &
                                rvanka%rvanka2DNavStOptC%p_KdiagonalM)
      else
        rvanka%rvanka2DNavStOptC%p_KdiagonalM => rvanka%rvanka2DNavStOptC%p_KldM
      end if

    end if

    if (lsysbl_isSubmatrixPresent (rmatrix,1,5)) then
      if (rmatrix%RmatrixBlock(1,5)%cmatrixFormat .eq. LSYSSC_MATRIXD) then
        call output_line ('Lumped mass matrices not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
        call sys_halt()
      end if

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,5),&
          rvanka%rvanka2DNavStOptC%p_DM15 )

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,4),&
          rvanka%rvanka2DNavStOptC%p_DM24 )

    end if

    ! Get the coupling matrix from the primal to the dual system.
    ! This may be the mass matrix or a completely decoupled matrix.
    ! In all cases, the submatrices have the same structure as the mass
    ! matrix, so if the structure of the mass matrix is not yet set,
    ! we have to fetch it!
    if (lsysbl_isSubmatrixPresent (rmatrix,4,1)) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,1),&
          rvanka%rvanka2DNavStOptC%p_DR41 )

      if (.not. associated(rvanka%rvanka2DNavStOptC%p_KcolM)) then
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(4,1),&
            rvanka%rvanka2DNavStOptC%p_KcolM)
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,1), &
            rvanka%rvanka2DNavStOptC%p_KldM )
      end if
    end if

    if (lsysbl_isSubmatrixPresent (rmatrix,4,2)) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,2),&
          rvanka%rvanka2DNavStOptC%p_DR42 )

      if (.not. associated(rvanka%rvanka2DNavStOptC%p_KcolM)) then
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(4,2),&
            rvanka%rvanka2DNavStOptC%p_KcolM)
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,2), &
            rvanka%rvanka2DNavStOptC%p_KldM )
      end if
    end if

    if (lsysbl_isSubmatrixPresent (rmatrix,5,1)) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(5,1),&
          rvanka%rvanka2DNavStOptC%p_DR51 )

      if (.not. associated(rvanka%rvanka2DNavStOptC%p_KcolM)) then
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(5,1),&
            rvanka%rvanka2DNavStOptC%p_KcolM)
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(5,1), &
            rvanka%rvanka2DNavStOptC%p_KldM )
      end if
    end if

    if (lsysbl_isSubmatrixPresent (rmatrix,5,2)) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(5,2),&
          rvanka%rvanka2DNavStOptC%p_DR52 )

      if (.not. associated(rvanka%rvanka2DNavStOptC%p_KcolM)) then
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(5,2),&
            rvanka%rvanka2DNavStOptC%p_KcolM)
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(5,2), &
            rvanka%rvanka2DNavStOptC%p_KldM )
      end if
    end if

    ! Get the multiplication factors of the submatrices.
    ! (-> for a later implementation; currently, the multipliers are not used!)
    rvanka%rvanka2DNavStOptC%Dmultipliers(1:6,1:6) = &
        rmatrix%RmatrixBlock(1:6,1:6)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest

    if (.not. associated(p_rblockDiscr)) then
      call output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    ! Get the discretisation structure of U,V and P from the block
    ! discretisation structure.
    ! We assume that the discretisation of the dual equations are the same
    ! as for the primal equations!
    rvanka%rvanka2DNavStOptC%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanka%rvanka2DNavStOptC%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanka%rvanka2DNavStOptC%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)

    if (rvanka%rvanka2DNavStOptC%p_rspatialDiscrU%inumFESpaces .ne. &
        rvanka%rvanka2DNavStOptC%p_rspatialDiscrV%inumFESpaces) then
      call output_line (&
          'Discretisation structures of X- and Y-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

    if ((rvanka%rvanka2DNavStOptC%p_rspatialDiscrP%inumFESpaces .ne. 1) .and. &
        (rvanka%rvanka2DNavStOptC%p_rspatialDiscrP%inumFESpaces .ne. &
          rvanka%rvanka2DNavStOptC%p_rspatialDiscrU%inumFESpaces)) then
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      call output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesOptC')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DNavierStokesOptC (rvanka2DNavStOptC, rvector, rrhs, domega, csubtype)

!<description>
  ! This routine applies the VANKA variant for 2D Navier-Stokes
  ! optimal control problems to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANKA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector

  ! The subtype of VANKA that should handle the above problem class.
  ! One of the VANKATP_xxxx constants, e.g. VANKATP_DIAGONAL.
  integer :: csubtype

!</input>

!<inputoutput>
  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavStOptC), intent(inout) :: rvanka2DNavStOptC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ielementdist
    integer, dimension(:), pointer :: p_IelementList
    type(t_elementDistribution), pointer :: p_relementDistrU
    type(t_elementDistribution), pointer :: p_relementDistrV
    type(t_elementDistribution), pointer :: p_relementDistrP

    ! 2D Navier Stokes problem.

    ! Loop through the element distributions of the velocity.
    do ielementdist = 1,rvanka2DNavStOptC%p_rspatialDiscrU%inumFESpaces

      ! Get the corresponding element distributions of U, V and P.
      p_relementDistrU => &
          rvanka2DNavStOptC%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanka2DNavStOptC%p_rspatialDiscrV%RelementDistr(ielementdist)

      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      if (rvanka2DNavStOptC%p_rspatialDiscrP%inumFESpaces .gt. 1) then
        p_relementDistrP => &
            rvanka2DNavStOptC%p_rspatialDiscrP%RelementDistr(ielementdist)
      else
        p_relementDistrP => &
            rvanka2DNavStOptC%p_rspatialDiscrP%RelementDistr(1)
      end if

      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y-velocity (and to that
      ! of the pressure).
      call storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)

      ! Which element combination do we have now?
      if ((elem_getPrimaryElement(p_relementDistrU%celement) .eq. EL_Q1T) .and. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .eq. EL_Q1T) .and. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .eq. EL_Q0)) then

        ! Q1~/Q1~/Q0 discretisation

        select case (csubtype)

        case (VANKATP_FULLOPTC_PRIMAL)

          ! Apply the conformal VANKA that allows different matrices
          ! in A11, A12, A21 and A22!
          call vanka_2DNSSOCQ1TQ0fullCoupCfFB (rvanka2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList,1)

        case (VANKATP_FULLOPTC_DUAL)

          ! Apply the conformal VANKA that allows different matrices
          ! in A11, A12, A21 and A22!
          call vanka_2DNSSOCQ1TQ0fullCoupCfFB (rvanka2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList,2)

        case (VANKATP_DIAGOPTC)

          ! Apply the conformal diagonal VANKA that allows different
          ! matrices in A11, A22, A33 and A44!
          call vanka_2DNSSOCQ1TQ0diagCoupConf (rvanka2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList)

        case default

          ! Apply the conformal VANKA that allows different matrices
          ! in A11, A12, A21 and A22!
          call vanka_2DNSSOCQ1TQ0fullCoupConf (rvanka2DNavStOptC, &
              rvector, rrhs, domega,p_IelementList)

        end select

      else

        call output_line ('Unsupported discretisation.',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNavierStokesOptC')
        call sys_halt()

      end if

    end do

  end subroutine

  ! ***************************************************************************
  ! 2D VANKA, 'diagonal' version for fully coupled Navier-Stokes systems with
  ! primal and dual equations.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A21  B1  M14          )
  !    ( A12  A22  B2       M25     )
  !    ( D1^T D2^T I1               )
  !    ( M41  M42      A44  A45  B1 )
  !    ( M51  M52      A54  A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! The mass matrices must all be the same except for their scaling factors!
  ! Here, a and b are arbitrary multiplication factors for the mass matrices.
  !
  ! I1 and I2 are two diagonal matrices in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! do not exist, what results in a '0' block in these positions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DNSSOCQ1TQ0diagCoupConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier Stokes optimal control problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DNSSOCQ1TQ0fullCoupConf supports fully coupled velocity submatrices.
  ! The matrices A11, A22, A44 and A55 must have the same structure.
  ! The matrices A12 and A21 must have the same structure.
  ! The matrices A45 and A54 must have the same structure.
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavStOptC), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof
    real(DP) :: dmult11,dmult22,dmult44,dmult55
    real(DP) :: dmultb1,dmultb2,dmultb3,dmultb4
    real(DP) :: dmultd1,dmultd2,dmultd3,dmultd4
    real(DP) :: dmult33,dmult66

    real(DP) :: dmult12,dmult21,dmult45,dmult54
    real(DP) :: dmult41,dmult52,dmult51,dmult42
    real(DP) :: dmult14,dmult25,dmult15,dmult24

    real(DP) :: di1,di2

    integer, dimension(:), pointer :: p_KcolA,p_KcolA12,p_KcolA45
    integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldA45,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA11,p_Da22,p_Da44,p_Da55
    real(DP), dimension(:), pointer :: p_DA12,p_Da21,p_Da45,p_Da54
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1,p_DB2
    real(DP), dimension(:), pointer :: p_DB3,p_DB4
    real(DP), dimension(:), pointer :: p_DD1,p_DD2
    real(DP), dimension(:), pointer :: p_DD3,p_DD4
    real(DP), dimension(:), pointer :: p_DA33,p_DA66
    real(DP), dimension(:), pointer :: p_DA14,p_Da25
    real(DP), dimension(:), pointer :: p_DA41,p_Da42,p_Da52,p_Da51
    real(DP), dimension(:), pointer :: p_DA15,p_Da24
    integer, dimension(:), pointer :: p_KdiagonalA33,p_KdiagonalA66

    logical :: bhaveA12,bhaveA45,bhaveA14,bhaveA41,bhaveA42,bhaveA24

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetu,ioffsetv,ioffsetp
    integer :: ioffsetl1,ioffsetl2,ioffsetxi
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 4
    integer, parameter :: lofsp = 8
    real(DP) :: daux1,daux2,daux3,daux4
    real(DP) :: FFp1,FFp2,FFd1,FFd2,FFpp,FFdp

    ! Local arrays for informations about one element -- for primal and dual space.
    real(DP), dimension(4) :: AA11,AA22,BB1,BB2,DD1,DD2
    real(DP), dimension(9) :: FFp,UUp
    real(DP), dimension(4) :: AA33,AA44,BB3,BB4,DD3,DD4
    real(DP), dimension(9) :: FFd,UUd
    integer, dimension(4) :: idofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA11
    p_KldA => rvanka%p_KldA11
    p_KdiagonalA => rvanka%p_KdiagonalA11
    p_KcolA12 => rvanka%p_KcolA12
    p_KldA12 => rvanka%p_KldA12
    p_KcolA45 => rvanka%p_KcolA45
    p_KldA45 => rvanka%p_KldA45
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22
    p_Da44 => rvanka%p_Da44
    p_Da55 => rvanka%p_Da55
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DB3 => rvanka%p_DB1
    p_DB4 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2
    p_DD3 => rvanka%p_DD1
    p_DD4 => rvanka%p_DD2
    p_DA33 => rvanka%p_DA33
    p_DA66 => rvanka%p_DA66
    p_KdiagonalA33 => rvanka%p_KdiagonalA33
    p_KdiagonalA66 => rvanka%p_KdiagonalA66

    p_Da12 => rvanka%p_Da12
    p_Da21 => rvanka%p_Da21
    p_Da45 => rvanka%p_Da45
    p_Da54 => rvanka%p_Da54
    p_DA14 => rvanka%p_Dm14
    p_DA15 => rvanka%p_Dm15
    p_DA24 => rvanka%p_Dm24
    p_Da25 => rvanka%p_Dm25
    p_DA41 => rvanka%p_Dr41
    p_Da42 => rvanka%p_Dr42
    p_Da52 => rvanka%p_Dr52
    p_Da51 => rvanka%p_Dr51

    dmult11 = rvanka%Dmultipliers(1,1)
    dmult22 = rvanka%Dmultipliers(2,2)
    dmult44 = rvanka%Dmultipliers(4,4)
    dmult55 = rvanka%Dmultipliers(5,5)

    dmultb1 = rvanka%Dmultipliers(1,3)
    dmultb2 = rvanka%Dmultipliers(2,3)
    dmultb3 = rvanka%Dmultipliers(4,6)
    dmultb4 = rvanka%Dmultipliers(5,6)

    dmultd1 = rvanka%Dmultipliers(3,1)
    dmultd2 = rvanka%Dmultipliers(3,2)
    dmultd3 = rvanka%Dmultipliers(6,4)
    dmultd4 = rvanka%Dmultipliers(6,5)

    dmult33 = rvanka%Dmultipliers(3,3)
    dmult66 = rvanka%Dmultipliers(6,6)

    dmult12 = rvanka%Dmultipliers(1,2)
    dmult21 = rvanka%Dmultipliers(2,1)
    dmult45 = rvanka%Dmultipliers(4,5)
    dmult54 = rvanka%Dmultipliers(5,4)
    dmult41 = rvanka%Dmultipliers(4,1)
    dmult52 = rvanka%Dmultipliers(5,2)
    dmult51 = rvanka%Dmultipliers(5,1)
    dmult42 = rvanka%Dmultipliers(4,2)
    dmult14 = rvanka%Dmultipliers(1,4)
    dmult25 = rvanka%Dmultipliers(2,5)
    dmult15 = rvanka%Dmultipliers(1,5)
    dmult24 = rvanka%Dmultipliers(2,4)

    bhaveA12 = associated(p_Da12) .and. (dmult12 .ne. 0.0_DP)
    bhaveA45 = associated(p_Da45) .and. (dmult45 .ne. 0.0_DP)
    bhaveA14 = associated(p_DA14) .and. (dmult14 .ne. 0.0_DP)
    bhaveA24 = associated(p_DA24) .and. (dmult24 .ne. 0.0_DP)
    bhaveA41 = associated(p_DA41) .and. (dmult41 .ne. 0.0_DP)
    bhaveA42 = associated(p_Da42) .and. (dmult42 .ne. 0.0_DP)

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the solution components
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetl1 = ioffsetp+rvector%RvectorBlock(3)%NEQ
    ioffsetl2 = ioffsetl1+rvector%RvectorBlock(4)%NEQ
    ioffsetxi = ioffsetl2+rvector%RvectorBlock(5)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP

      FFpp = p_Drhs(iel+ioffsetp)
      FFdp = p_Drhs(iel+ioffsetxi)

      ! Loop over all 4 U-nodes of that element.
      do inode=1,4

        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA11(inode) = dmult11*p_DA11(p_KdiagonalA(idof))
        AA22(inode) = dmult22*p_DA22(p_KdiagonalA(idof))
        AA33(inode) = dmult44*p_Da44(p_KdiagonalA(idof))
        AA44(inode) = dmult55*p_Da55(p_KdiagonalA(idof))


        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FFp1 = p_Drhs(idof+ioffsetu)
        FFp2 = p_Drhs(idof+ioffsetv)

        FFd1 = p_Drhs(idof+ioffsetl1)
        FFd2 = p_Drhs(idof+ioffsetl2)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the four velocity unknowns and the
        ! pressure unknown on the current element; these five unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity
        ! components and B being an (2*4) x 1 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |        |             |        |        |
        !   --|---X----|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        daux3 = 0.0_DP
        daux4 = 0.0_DP
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux1 = daux1 + p_Da11(ia)*p_Dvector(J+ioffsetu)
          daux2 = daux2 + p_Da22(ia)*p_Dvector(J+ioffsetv)
          daux3 = daux3 + p_Da44(ia)*p_Dvector(J+ioffsetl1)
          daux4 = daux4 + p_Da55(ia)*p_Dvector(J+ioffsetl2)
        end do
        FFp1 = FFp1 - dmult11*daux1
        FFp2 = FFp2 - dmult22*daux2
        FFd1 = FFd1 - dmult44*daux3
        FFd2 = FFd2 - dmult55*daux4

        ! There are probably some more defects to calculate.
        if (bhaveA14) then
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do ia = ia1,ia2
            J = p_KcolA(ia)
            daux1 = daux1 + p_DA14(ia)*p_Dvector(J+ioffsetl1)
            daux2 = daux2 + p_DA25(ia)*p_Dvector(J+ioffsetl2)
          end do
          FFp1 = FFp1 - dmult14*daux1
          FFp2 = FFp2 - dmult25*daux2
        end if

        if (bhaveA41) then
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do ia = ia1,ia2
            J = p_KcolA(ia)
            daux1 = daux1 + p_DA41(ia)*p_Dvector(J+ioffsetu)
            daux2 = daux2 + p_DA52(ia)*p_Dvector(J+ioffsetv)
          end do
          FFd1 = FFd1 - dmult41*daux1
          FFd2 = FFd2 - dmult52*daux2
        end if

        if (bhaveA12) then
          ia1 = p_KldA12(idof)
          ia2 = p_KldA12(idof+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do ia = ia1,ia2
            J = p_KcolA(ia)
            daux1 = daux1 + p_DA12(ia)*p_Dvector(J+ioffsetv)
            daux2 = daux2 + p_DA21(ia)*p_Dvector(J+ioffsetu)
          end do
          FFp1 = FFp1 - dmult12*daux1
          FFp2 = FFp2 - dmult21*daux2
        end if

        if (bhaveA45) then
          ia1 = p_KldA45(idof)
          ia2 = p_KldA45(idof+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do ia = ia1,ia2
            J = p_KcolA(ia)
            daux1 = daux1 + p_DA45(ia)*p_Dvector(J+ioffsetl2)
            daux2 = daux2 + p_DA54(ia)*p_Dvector(J+ioffsetl1)
          end do
          FFd1 = FFd1 - dmult45*daux1
          FFd2 = FFd2 - dmult54*daux2
        end if

        if (bhaveA42) then
          ia1 = p_KldA45(idof)
          ia2 = p_KldA45(idof+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do ia = ia1,ia2
            J = p_KcolA(ia)
            daux1 = daux1 + p_DA42(ia)*p_Dvector(J+ioffsetv)
            daux2 = daux2 + p_DA51(ia)*p_Dvector(J+ioffsetu)
          end do
          FFd1 = FFd1 - dmult42*daux1
          FFd2 = FFd2 - dmult51*daux2
        end if

        if (bhaveA24) then
          ia1 = p_KldA45(idof)
          ia2 = p_KldA45(idof+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do ia = ia1,ia2
            J = p_KcolA(ia)
            daux1 = daux1 + p_DA15(ia)*p_Dvector(J+ioffsetl2)
            daux2 = daux2 + p_DA24(ia)*p_Dvector(J+ioffsetl1)
          end do
          FFp1 = FFp1 - dmult15*daux1
          FFp2 = FFp2 - dmult24*daux2
        end if

        ! Finally subtract B*p: f_i = (f_i-Aui) - Bi pi
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux1 = p_Dvector(j+ioffsetp)
          daux2 = p_Dvector(j+ioffsetxi)

          FFp1 = FFp1-dmultb1*p_DB1(ib)*daux1
          FFp2 = FFp2-dmultb2*p_DB2(ib)*daux1

          FFd1 = FFd1-dmultb3*p_DB3(ib)*daux2
          FFd2 = FFd2-dmultb4*p_DB4(ib)*daux2
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = dmultb1*p_DB1(ib)
            BB2(inode) = dmultb2*p_DB2(ib)

            BB3(inode) = dmultb3*p_DB3(ib)
            BB4(inode) = dmultb4*p_DB4(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = dmultd1*p_DD1(ib)
            DD2(inode) = dmultd2*p_DD2(ib)

            DD3(inode) = dmultd3*p_DD3(ib)
            DD4(inode) = dmultd4*p_DD4(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FFpp = FFpp &
                 - DD1(inode)*p_Dvector(idof+ioffsetu) &
                 - DD2(inode)*p_Dvector(idof+ioffsetv)

            FFdp = FFdp &
                 - DD3(inode)*p_Dvector(idof+ioffsetl1) &
                 - DD4(inode)*p_Dvector(idof+ioffsetl2)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

        FFp(inode) = FFp1
        FFp(inode+lofsv) = FFp2
        FFd(inode) = FFd1
        FFd(inode+lofsv) = FFd2
        FFp(1+lofsp) = FFpp
        FFd(1+lofsp) = FFdp

      end do ! inode

      ! If we have blocks at A33 or A66, get the corresponding elements.
      ! We need the IF-commands here as the arrays p_DaXX/p_KdiagonalXX may
      ! be undefined.

      di1 = 0.0_DP
      di2 = 0.0_DP

      if (dmult33 .ne. 0.0_DP) then
        di1 = dmult33*p_DA33(p_KdiagonalA33(iel))
      end if

      if (dmult66 .ne. 0.0_DP) then
        di2 = dmult66*p_DA66(p_KdiagonalA66(iel))
      end if

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_2DSPQ1TQ0simple2 (UUp,FFp,AA11,AA22,BB1,BB2,DD1,DD2,di1)
      call vanka_getcorr_2DSPQ1TQ0simple2 (UUd,FFd,AA33,AA44,BB3,BB4,DD3,DD4,di2)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,4
        p_Dvector(idofGlobal(inode)+ioffsetu) &
          = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * UUp(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UUp(inode+lofsv)

        p_Dvector(idofGlobal(inode)+ioffsetl1) &
          = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * UUd(inode)
        p_Dvector(idofGlobal(inode)+ioffsetl2) &
          = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * UUd(inode+lofsv)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UUp(1+lofsp)

      p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + domega * UUd(1+lofsp)

    end do ! iel

  end subroutine

  ! ***************************************************************************
  ! 2D VANKA, 'full' version for fully coupled Navier-Stokes systems with
  ! primal and dual equations.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1  aM           )
  !    ( A21  A22  B2       aM      )
  !    ( D1^T D2^T I1               )
  !    ( bR            A44  A45  B1 )
  !    (      bR       A54  A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! The mass matrices must all be the same except for their scaling factors!
  ! Here, a and b are arbitrary multiplication factors for the mass matrices.
  !
  ! I1 and I2 are two diagonal matrices in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! do not exist, what results in a '0' block in these positions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DNSSOCQ1TQ0fullCoupConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier Stokes optimal control problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DNSSOCQ1TQ0fullCoupConf supports fully coupled velocity submatrices.
  ! The matrices A11, A22, A44 and A55 must have the same structure.
  ! The matrices A12 and A21 must have the same structure.
  ! The matrices A45 and A54 must have the same structure.
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavStOptC), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA11
    integer, dimension(:), pointer :: p_KldA11
    real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolA45
    integer, dimension(:), pointer :: p_KldA45
    real(DP), dimension(:), pointer :: p_DA44,p_DA45,p_DA54,p_DA55
    real(DP), dimension(:), pointer :: p_DR41,p_DR52,p_DR51,p_DR42
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12
    integer, dimension(:), pointer :: p_KcolM
    integer, dimension(:), pointer :: p_KldM
    real(DP), dimension(:), pointer :: p_DM
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2
    real(DP), dimension(:), pointer :: p_Da33,p_Da66
    integer, dimension(:), pointer :: p_KdiagonalA33,p_KdiagonalA66

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 4      ! Q1T = 4 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! QQ0 = 1 DOF`s per pressure
    integer, parameter :: nndualvel = 4      ! Q1T = 4 DOF`s per dual velocity
    integer, parameter :: nndualpressure = 1 ! QQ0 = 1 DOF`s per dual pressure
    integer, parameter :: nnprimal = 2*nnvel+nnpressure ! Q1~/Q1~/Q0 = 4+4+1 = 9 DOF`s per element
    integer, parameter :: nnld = 2*nnprimal
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! Offsets of the 'local' solution parts in the 'local' solution vector
    integer, parameter :: lofsu = 0
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    integer, parameter :: lofsl1 = 2*nnvel+1
    integer, parameter :: lofsl2 = 2*nnvel+1+nnvel
    integer, parameter :: lofsxi = 2*nnvel+1+2*nnvel

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! Offset information in arrays.
    ! Primal variables
    integer :: ioffsetu,ioffsetv,ioffsetp,j

    ! Dual variables
    integer :: ioffsetl1,ioffsetl2,ioffsetxi

    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    real(DP) :: daux,daux2

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA11 => rvanka%p_KcolA11
    p_KldA11 => rvanka%p_KldA11
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22

    ! Structure of A12 is assumed to be the same as A21.
    ! Get A12 and A21 -- except for if the multipliers are =0, then
    ! we switch them off by nullifying the pointers.
    if (rvanka%Dmultipliers(1,2) .ne. 0.0_DP) then
      p_KcolA12 => rvanka%p_KcolA12
      p_KldA12 => rvanka%p_KldA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    else
      nullify(p_KcolA12)
      nullify(p_KldA12)
      nullify(p_DA12 )
      nullify(p_DA21 )
    end if

    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Structure of A44 is assumed to be the same as A55, A11 and A22
    p_DA44 => rvanka%p_DA44
    p_DA55 => rvanka%p_DA55

    ! Structure of A45 is assumed to be the same as A54
    p_KcolA45 => rvanka%p_KcolA45
    p_KldA45 => rvanka%p_KldA45
    p_DA45 => rvanka%p_DA45
    p_DA54 => rvanka%p_DA54

    ! Mass matrix - if it is given, otherwise the pointers will be set to NULL
    ! because of the initialisation of the structure!
    p_KcolM => rvanka%p_KcolM
    p_KldM => rvanka%p_KldM
    p_DM => rvanka%p_DM14

    ! Coupling matrix in the dual equation at position (4:5,1:2). For a standard
    ! system, there is A(4,1) = A(5,2) = M and A(5,1) = A(4,2) = 0.
    ! For a Newton system, this block is completely decoupled!
    p_DR41 => rvanka%p_DR41
    p_DR42 => rvanka%p_DR42
    p_DR51 => rvanka%p_DR51
    p_DR52 => rvanka%p_DR52

    ! Diagonal submatrices A33 and A66 (if they exist)
    if (rvanka%Dmultipliers(3,3) .ne. 0.0_DP) then
      p_Da33 => rvanka%p_DA33
      p_KdiagonalA33 => rvanka%p_KdiagonalA33
    else
      nullify(p_Da33)
      nullify(p_KdiagonalA33)
    end if

    if (rvanka%Dmultipliers(6,6) .ne. 0.0_DP) then
      p_Da66 => rvanka%p_DA66
      p_KdiagonalA66 => rvanka%p_KdiagonalA66
    else
      nullify(p_Da66)
      nullify(p_KdiagonalA66)
    end if

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    ! Get the offsets of lambda1, lambda2 and xi, so the offsets
    ! of the dual solution vectors.
    ioffsetl1 = ioffsetp+rvector%RvectorBlock(3)%NEQ
    ioffsetl2 = ioffsetl1+rvector%RvectorBlock(4)%NEQ
    ioffsetxi = ioffsetl2+rvector%RvectorBlock(5)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          4    P    2
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !
      !
      ! Fetch the pressure P on the current element into FF.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      ! Get the primal pressure
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the dual pressure
      FF(1+lofsxi) = p_Drhs(iel+ioffsetxi)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..4 = edge-NVT.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode.

        ! Primal equation
        FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! dual equation
        FF(inode+lofsl1) = p_Drhs(idof+ioffsetl1)
        FF(inode+lofsl2) = p_Drhs(idof+ioffsetl2)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B  M     ] (u)  = (f  )
        !    [ B^t 0        ] (p)    (g  )
        !    [ M      A   B ] (l)    (fl )
        !    [        B^t 0 ] (xi)   (fxi)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in rectangular
        ! systems of the form
        !
        !    [ === A^ === B~  === M^ ====   ] (| ) = (f1 )
        !    [ B~^t       I1~               ] (u )   (f2 )
        !    [ === M^ ===     === A^ === B~ ] (| )   (g  )
        !    [                B~^t       I2~] (p )   (fl1)
        !                                     (| )   (fl2)
        !                                     (l )   (flg)
        !                                     (| )
        !                                     (xi)
        !
        !
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! 2*1 pressure elements on the adjacent cells, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL.
        !
        !  (d1 ) = (f1 ) -  [ === A^ === B~  === M^ ====   ] (| )
        !  (d2 )   (f2 )    [ B~^t       I1~               ] (u )
        !  (dg )   (g  )    [ === M^ ===     === A^ === B~ ] (| )
        !  (dl1)   (fl1)    [                B~^t       I2~] (p )
        !  (dl2)   (fl2)                                     (| )
        !  (dlg)   (flg)                                     (l )
        !                                                    (| )
        !                                                    (xi)
        !
        ! Extract those entries in the A-, B- and M-matrices to our local
        ! matrix AA, which belong to DOF`s in our current solution vector.
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA11(idof)
        ia2 = p_KldA11(idof+1)-1
        do ia = ia1,ia2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - ( A11  .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .  A22  .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .  A44  .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .  A55  .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

          J = p_KcolA11(ia)

          ! Primal equation:
          FF(inode+lofsu) = FF(inode+lofsu) &
                          - rvanka%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
          FF(inode+lofsv) = FF(inode+lofsv) &
                          - rvanka%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! dual equation
          FF(inode+lofsl1) = FF(inode+lofsl1) &
                           - rvanka%Dmultipliers(4,4)*p_DA44(ia)*p_Dvector(J+ioffsetl1)
          FF(inode+lofsl2) = FF(inode+lofsl2) &
                           - rvanka%Dmultipliers(5,5)*p_DA55(ia)*p_Dvector(J+ioffsetl2)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanka%Dmultipliers(1,1)
              AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanka%Dmultipliers(2,2)

              AA (inode+lofsl1,k+lofsl1) = p_DA44(ia)*rvanka%Dmultipliers(4,4)
              AA (inode+lofsl2,k+lofsl2) = p_DA55(ia)*rvanka%Dmultipliers(5,5)
              exit
            end if
          end do
        end do

        ! Process the 'off-diagonal' matrices A12 and A21

        if (associated(p_KldA12)) then
          ia1 = p_KldA12(idof)
          ia2 = p_KldA12(idof+1)-1
          do ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .  A12  .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   ( A21  .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA12(ia)
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanka%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanka%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanka%Dmultipliers(1,2)
                AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanka%Dmultipliers(2,1)
                exit
              end if
            end do
          end do
        end if

        ! Process the 'off-diagonal' matrices A45 and A54 if they exist
        if (associated(p_KldA45)) then
          ia1 = p_KldA45(idof)
          ia2 = p_KldA45(idof+1)-1
          do ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .  A45  .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .  A54  .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA45(ia)
            FF(inode+lofsl1) = FF(inode+lofsl1) &
                             - rvanka%Dmultipliers(4,5)*p_DA45(ia)*p_Dvector(J+ioffsetl2)
            FF(inode+lofsl2) = FF(inode+lofsl2) &
                             - rvanka%Dmultipliers(5,4)*p_DA54(ia)*p_Dvector(J+ioffsetl1)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsl1,k+lofsl2) = p_DA45(ia)*rvanka%Dmultipliers(4,5)
                AA (inode+lofsl2,k+lofsl1) = p_DA54(ia)*rvanka%Dmultipliers(5,4)
                exit
              end if
            end do
          end do
        end if

        ! Process A33 if it exists

        if (associated(p_KdiagonalA33)) then

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   I1  .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
          !
          ! IEL is the pressure DOF which we have to tackle.

          daux = rvanka%Dmultipliers(3,3)
          FF(1+lofsp) = FF(1+lofsp) &
                      - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
          AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

        end if

        ! Process A66 if it exists

        if (associated(p_KdiagonalA66)) then

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   I2 ) ( xi )
          !
          ! IEL is the pressure DOF which we have to tackle.

          daux = rvanka%Dmultipliers(6,6)
          FF(1+lofsxi) = FF(1+lofsxi) &
                       - daux*p_DA66(p_KdiagonalA66(IEL))*p_Dvector(IEL+ioffsetxi)
          AA(1+lofsxi,1+lofsxi) = daux*p_DA66(p_KdiagonalA66(IEL))

        end if

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .  B1   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .  B2   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .  B1  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .  B2  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

          J = p_KcolB(ib)

          ! primal equation
          daux = p_Dvector(j+ioffsetp)
          FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanka%Dmultipliers(1,3)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanka%Dmultipliers(2,3)

          ! dual equation
          daux2 = p_Dvector(j+ioffsetxi)
          FF(inode+lofsl1) = FF(inode+lofsl1)-p_DB1(ib)*daux2 * rvanka%Dmultipliers(4,6)
          FF(inode+lofsl2) = FF(inode+lofsl2)-p_DB2(ib)*daux2 * rvanka%Dmultipliers(5,6)

          ! Do not incorporate the B-matrices into AA yet; this will come later!
        end do

        ! The mass matrix defect.
        if (associated(p_DM)) then
          ! We assume: multiplier of A(1,4) = multiplier of A(2,5)
          daux = rvanka%Dmultipliers(1,4)

          ! We assume: multiplier of A(4,1) = multiplier of A(5,2)
          daux2 = rvanka%Dmultipliers(4,1)

          ia1 = p_KldM(idof)
          ia2 = p_KldM(idof+1)-1
          do ia = ia1,ia2

            J = p_KcolM(ia)

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .  aM   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .  aM   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            FF(inode+lofsu) = FF(inode+lofsu)-daux*p_DM(ia)*p_Dvector(J+ioffsetl1)
            FF(inode+lofsv) = FF(inode+lofsv)-daux*p_DM(ia)*p_Dvector(J+ioffsetl2)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsu,k+lofsl1) = daux*p_DM(ia)
                AA (inode+lofsv,k+lofsl2) = daux*p_DM(ia)
                exit
              end if
            end do
          end do
        end if

        ! The defect in the coupling matrix from the primal to the dual system
        if (associated(p_DR41)) then
          ! Get the multipliers
          daux = rvanka%Dmultipliers(4,1)
          daux2 = rvanka%Dmultipliers(5,2)

          ia1 = p_KldM(idof)
          ia2 = p_KldM(idof+1)-1
          do ia = ia1,ia2

            J = p_KcolM(ia)

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   ( bR   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .  bR   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            FF(inode+lofsl1) = FF(inode+lofsl1)-daux*p_DR41(ia)*p_Dvector(J+ioffsetu)
            FF(inode+lofsl2) = FF(inode+lofsl2)-daux2*p_DR52(ia)*p_Dvector(J+ioffsetv)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsl1,k+lofsu) = daux*p_DR41(ia)
                AA (inode+lofsl2,k+lofsv) = daux2*p_DR52(ia)
                exit
              end if
            end do
          end do
        end if

        if (associated(p_DR51)) then
          ! Get the multipliers
          daux = rvanka%Dmultipliers(5,1)
          daux2 = rvanka%Dmultipliers(4,2)

          ia1 = p_KldM(idof)
          ia2 = p_KldM(idof+1)-1
          do ia = ia1,ia2

            J = p_KcolM(ia)

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .  bR   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   ( bR   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            FF(inode+lofsl1) = FF(inode+lofsl1)-daux2*p_DR42(ia)*p_Dvector(J+ioffsetv)
            FF(inode+lofsl2) = FF(inode+lofsl2)-daux*p_DR51(ia)*p_Dvector(J+ioffsetu)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsl2,k+lofsu) = daux*p_DR51(ia)
                AA (inode+lofsl1,k+lofsv) = daux2*p_DR42(ia)
                exit
              end if
            end do
          end do
        end if

        ! THe next loop will determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
          !   ( dp  )   ( dp  )   ( D1  D2   .   .   .   .  ) ( p  )
          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
          !   ( dxi )   ( dxi )   (  .   .   .  D1  D2   .  ) ( xi )
          !
          ! In AA, we simultaneously set up (locally):
          !
          !   (  .   .  B1   .   .   .  )
          !   (  .   .  B2   .   .   .  )
          !   ( D1  D2   .   .   .   .  )
          !   (  .   .   .   .   .  B1  )
          !   (  .   .   .   .   .  B2  )
          !   (  .   .   .  D1  D2   .  )

          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices.
            ! Primal equation
            AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanka%Dmultipliers(1,3)
            AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanka%Dmultipliers(2,3)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanka%Dmultipliers(3,1)
            AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanka%Dmultipliers(3,2)

            ! The same for the dual equation
            AA(inode+lofsl1,1+lofsxi) = p_DB1(ib) * rvanka%Dmultipliers(4,6)
            AA(inode+lofsl2,1+lofsxi) = p_DB2(ib) * rvanka%Dmultipliers(5,6)

            AA(1+lofsxi,inode+lofsl1) = p_DD1(ib) * rvanka%Dmultipliers(6,4)
            AA(1+lofsxi,inode+lofsl2) = p_DD2(ib) * rvanka%Dmultipliers(6,5)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                        - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)

            ! The same for the dual pressure
            FF(1+lofsxi) = FF(1+lofsxi) &
                         - AA(1+lofsxi,inode+lofsl1)*p_Dvector(idof+ioffsetl1) &
                         - AA(1+lofsxi,inode+lofsl2)*p_Dvector(idof+ioffsetl2)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! For C, we use our local AA, i.e. applying C^{-1} means to
      ! solve the local system AA dd = FF for dd. The local defect dd is then
      ! added back to the global solution vector.

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y

        do inode=1,nnvel
          ! Update of the primal velocity vectors
          p_Dvector(idofGlobal(inode)+ioffsetu) &
            = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)

          ! Update of the dual velocity vectors
          p_Dvector(idofGlobal(inode)+ioffsetl1) &
            = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * FF(inode+lofsl1)
          p_Dvector(idofGlobal(inode)+ioffsetl2) &
            = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * FF(inode+lofsl2)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)

        p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + &
                                  domega * FF(1+lofsxi)

      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNSSOCQ1TQ0fullCoupConf')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************
  ! 2D VANKA, 'full' version for fully coupled Navier-Stokes systems with
  ! primal and dual equations. Forward-Backward-variant that handles only either
  ! the primal or the dual system.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    ( A11  A12  B1  aM           )
  !    ( A21  A22  B2       aM      )
  !    ( D1^T D2^T I1               )
  !    ( bM            A44  A45  B1 )
  !    (      bM       A54  A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! The mass matrices must all be the same except for their scaling factors!
  ! Here, a and b are arbitrary multiplication factors for the mass matrices.
  !
  ! I1 and I2 are two diagonal matrices in format 9, which may or may not
  ! exist in the system. For usual saddle point problems, these matrices
  ! do not exist, what results in a '0' block in these positions.
  !
  ! This variant decouples the primal system from the dual one. Depending on
  ! the parameter csystemType in the VANKA structure, only parts of the system
  ! are processed:
  !
  ! csystemType = 1: Apply VANKA to the system
  !
  !    ( A11  A12  B1  aM           )
  !    ( A21  A22  B2       aM      )
  !    ( D1^T D2^T I1               )
  !    (               I            )
  !    (                    I       )
  !    (                         I  )
  !
  ! csystemType = 2: Apply VANKA to the system
  !
  !    ( I                          )
  !    (      I                     )
  !    (           I                )
  !    ( bM            A44  A45  B1 )
  !    (      bM       A54  A55  B2 )
  !    (               D1^T D2^T I2 )
  !
  ! This variant is typically used in a nonstationary environment as Gauss-Seidel
  ! preconditioner which is at first applied to the primal equation before being
  ! applied to the dual equation.
  !
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DNSSOCQ1TQ0fullCoupCfFB (rvanka, rvector, rrhs, domega, IelementList,&
      csystemPart)

!<description>
  ! This routine applies the specialised forward-backward full local system
  ! VANKA algorithm for 2D Navier Stokes optimal control problems with
  ! Q1~/Q0 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DNSSOCQ1TQ0fullCoupCfFB supports fully coupled velocity submatrices.
  ! The matrices A11, A22, A44 and A55 must have the same structure.
  ! The matrices A12 and A21 must have the same structure.
  ! The matrices A45 and A54 must have the same structure.
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavStOptC), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList

  ! Identifier for the part of the equation, where VANKA should be
  ! applied.
  ! =0: apply VANKA to the whole primal-dual-coupled system
  ! =1: apply VANKA only to the primal system, take the dual one as constant
  ! =2: apply VANKA only to the dual system, take the primal one as constant
  integer, intent(in) :: csystemPart

!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA11
    integer, dimension(:), pointer :: p_KldA11
    real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolA45
    integer, dimension(:), pointer :: p_KldA45
    real(DP), dimension(:), pointer :: p_DA44,p_DA45,p_DA54,p_DA55
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12
    integer, dimension(:), pointer :: p_KcolM
    integer, dimension(:), pointer :: p_KldM
    real(DP), dimension(:), pointer :: p_DM
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2
    real(DP), dimension(:), pointer :: p_Da33,p_Da66
    integer, dimension(:), pointer :: p_KdiagonalA33,p_KdiagonalA66

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 4      ! Q1T = 4 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! QQ0 = 1 DOF`s per pressure
    integer, parameter :: nndualvel = 4      ! Q1T = 4 DOF`s per dual velocity
    integer, parameter :: nndualpressure = 1 ! QQ0 = 1 DOF`s per dual pressure
    integer, parameter :: nnprimal = 2*nnvel+nnpressure ! Q1~/Q1~/Q0 = 4+4+1 = 9 DOF`s per element
    integer, parameter :: nnld = 2*nnprimal
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! Offsets of the 'local' solution parts in the 'local' solution vector
    integer, parameter :: lofsu = 0
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    integer, parameter :: lofsl1 = 2*nnvel+1
    integer, parameter :: lofsl2 = 2*nnvel+1+nnvel
    integer, parameter :: lofsxi = 2*nnvel+1+2*nnvel

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! Offset information in arrays.
    ! Primal variables
    integer :: ioffsetu,ioffsetv,ioffsetp,j

    ! Dual variables
    integer :: ioffsetl1,ioffsetl2,ioffsetxi

    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    real(DP) :: daux,daux2

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA11 => rvanka%p_KcolA11
    p_KldA11 => rvanka%p_KldA11
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22

    ! Structure of A12 is assumed to be the same as A21.
    ! Get A12 and A21 -- except for if the multipliers are =0, then
    ! we switch them off by nullifying the pointers.
    if (rvanka%Dmultipliers(1,2) .ne. 0.0_DP) then
      p_KcolA12 => rvanka%p_KcolA12
      p_KldA12 => rvanka%p_KldA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    else
      nullify(p_KcolA12)
      nullify(p_KldA12)
      nullify(p_DA12 )
      nullify(p_DA21 )
    end if

    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Structure of A44 is assumed to be the same as A55, A11 and A22
    p_DA44 => rvanka%p_DA44
    p_DA55 => rvanka%p_DA55

    ! Structure of A45 is assumed to be the same as A54
    p_KcolA45 => rvanka%p_KcolA45
    p_KldA45 => rvanka%p_KldA45
    p_DA45 => rvanka%p_DA45
    p_DA54 => rvanka%p_DA54

    ! Mass matrix - if it is given, otherwise the pointers will be set to NULL
    ! because of the initialisation of the structure!
    p_KcolM => rvanka%p_KcolM
    p_KldM => rvanka%p_KldM
    p_DM => rvanka%p_DM14

    ! Diagonal submatrices A33 and A66 (if they exist)
    if (rvanka%Dmultipliers(3,3) .ne. 0.0_DP) then
      p_Da33 => rvanka%p_DA33
      p_KdiagonalA33 => rvanka%p_KdiagonalA33
    else
      nullify(p_Da33)
      nullify(p_KdiagonalA33)
    end if

    if (rvanka%Dmultipliers(6,6) .ne. 0.0_DP) then
      p_Da66 => rvanka%p_DA66
      p_KdiagonalA66 => rvanka%p_KdiagonalA66
    else
      nullify(p_Da66)
      nullify(p_KdiagonalA66)
    end if

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    ! Get the offsets of lambda1, lambda2 and xi, so the offsets
    ! of the dual solution vectors.
    ioffsetl1 = ioffsetp+rvector%RvectorBlock(3)%NEQ
    ioffsetl2 = ioffsetl1+rvector%RvectorBlock(4)%NEQ
    ioffsetxi = ioffsetl2+rvector%RvectorBlock(5)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    select case (csystemPart)
    ! -------------------------------------------------------------------------
    case (0)
      ! VANKA applied to full 18x18 system

      do ielidx=1,size(IelementList)

        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)

        ! Clear the 'local system matrix'.
        AA(:,:) = 0.0_DP

        ! We now have the element
        !
        ! +---------+                       +----3----+
        ! |         |                       |         |
        ! |   IEL   |   with DOF`s          4    P    2
        ! |         |                       |    Q0   |
        ! +---------+                       +----1----+
        !
        !
        ! Fetch the pressure P on the current element into FF.
        ! The numbers of the DOF`s coincide with the definition
        ! in dofmapping.f90!

        ! Get the primal pressure
        FF(1+lofsp) = p_Drhs(iel+ioffsetp)

        ! Get the dual pressure
        FF(1+lofsxi) = p_Drhs(iel+ioffsetxi)

        ! Get the velocity DOF`s on the current element.
        ! We assume: DOF 1..4 = edge.
        ! That is the same implementation as in dofmapping.f90!
        IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

        ! Loop over all U-nodes of that element.
        do inode=1,nnvel

          ! Get the DOF we have to tackle:
          idof = IdofGlobal(inode)

          ! Set FF initially to the value of the right hand
          ! side vector that belongs to our current DOF corresponding
          ! to inode.

          ! Primal equation
          FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
          FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

          ! dual equation
          FF(inode+lofsl1) = p_Drhs(idof+ioffsetl1)
          FF(inode+lofsl2) = p_Drhs(idof+ioffsetl2)

          ! What do we have at this point?
          ! FF     : "local" RHS vector belonging to the DOF`s on the
          !          current element
          ! AA     : Diagonal entries of A belonging to these DOF`s
          !
          ! And at the moment:
          ! idof      : number of current DOF on element IEL
          ! inode     : "local" number of DOF on element IEL, i.e.
          !              number of the edge
          !
          ! Now comes the crucial point with the "update": How to
          ! subsequently update the vertex values, such that the whole
          ! thing still converges to the solution, even if a node
          ! is updated more than once? Here, we use a typical
          ! matrix-decomposition approach:
          !
          ! Again consider the problem:
          !
          !    [ A   B  M     ] (u)  = (f  )
          !    [ B^t 0        ] (p)    (g  )
          !    [ M      A   B ] (l)    (fl )
          !    [        B^t 0 ] (xi)   (fxi)
          !
          ! We assume, that all components in the vector (u,p) are
          ! given - except for the velocity and pressure unknowns
          ! on the current element; these 21 unknowns
          ! are located anywhere in the (u,p) vector. The idea is to
          ! shift "everything known" to the right hand side to obtain
          ! a system for only these unknowns!
          !
          ! Extracting all the lines of the system that correspond to
          ! DOF`s on our single element IEL results in rectangular
          ! systems of the form
          !
          !    [ === A^ === B~  === M^ ====   ] (| ) = (f1 )
          !    [ B~^t       I1~               ] (u )   (f2 )
          !    [ === M^ ===     === A^ === B~ ] (| )   (g  )
          !    [                B~^t       I2~] (p )   (fl1)
          !                                     (| )   (fl2)
          !                                     (l )   (flg)
          !                                     (| )
          !                                     (xi)
          !
          !
          ! B~ is a 8 x 2 matrix: As every velocity couples with at most
          ! 2*1 pressure elements on the adjacent cells, so we have
          ! 2 columns in the B-matrix.
          !
          !        IEL                              IEL
          !     |--------|             |--------|--------|
          !     |        |             |        |        |
          !     |   P    |      or     |   Q    X   P    |
          !     |   X    |             |        |        |
          !   --|--------|--           |--------|--------|
          !
          !
          ! Now, throw all summands to the RHS vector to build a local
          ! 'defect' on our single element IEL.
          !
          !  (d1 ) = (f1 ) -  [ === A^ === B~  === M^ ====   ] (| )
          !  (d2 )   (f2 )    [ B~^t       I1~               ] (u )
          !  (dg )   (g  )    [ === M^ ===     === A^ === B~ ] (| )
          !  (dl1)   (fl1)    [                B~^t       I2~] (p )
          !  (dl2)   (fl2)                                     (| )
          !  (dlg)   (flg)                                     (l )
          !                                                    (| )
          !                                                    (xi)
          !
          ! Extract those entries in the A-, B- and M-matrices to our local
          ! matrix AA, which belong to DOF`s in our current solution vector.
          !
          ! At first build: fi = fi-Aui

          ia1 = p_KldA11(idof)
          ia2 = p_KldA11(idof+1)-1
          do ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - ( A11  .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .  A22  .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .  A44  .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .  A55  .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA11(ia)

            ! Primal equation:
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanka%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanka%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

            ! dual equation
            FF(inode+lofsl1) = FF(inode+lofsl1) &
                            - rvanka%Dmultipliers(4,4)*p_DA44(ia)*p_Dvector(J+ioffsetl1)
            FF(inode+lofsl2) = FF(inode+lofsl2) &
                            - rvanka%Dmultipliers(5,5)*p_DA55(ia)*p_Dvector(J+ioffsetl2)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanka%Dmultipliers(1,1)
                AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanka%Dmultipliers(2,2)

                AA (inode+lofsl1,k+lofsl1) = p_DA44(ia)*rvanka%Dmultipliers(4,4)
                AA (inode+lofsl2,k+lofsl2) = p_DA55(ia)*rvanka%Dmultipliers(5,5)
                exit
              end if
            end do
          end do

          ! Process the 'off-diagonal' matrices A12 and A21

          if (associated(p_KldA12)) then
            ia1 = p_KldA12(idof)
            ia2 = p_KldA12(idof+1)-1
            do ia = ia1,ia2
              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .  A12  .   .   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   ( A21  .   .   .   .   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              J = p_KcolA12(ia)
              FF(inode+lofsu) = FF(inode+lofsu) &
                              - rvanka%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
              FF(inode+lofsv) = FF(inode+lofsv) &
                              - rvanka%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)

              ! Whereever we find a DOF that couples to another DOF on the
              ! same element, we put that to both A-blocks of our local matrix.
              do k=1,nnvel
                if (j .eq. IdofGlobal(k)) then
                  AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanka%Dmultipliers(1,2)
                  AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanka%Dmultipliers(2,1)
                  exit
                end if
              end do
            end do
          end if

          ! Process the 'off-diagonal' matrices A45 and A54 if they exist
          if (associated(p_KldA45)) then
            ia1 = p_KldA45(idof)
            ia2 = p_KldA45(idof+1)-1
            do ia = ia1,ia2
              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !   ( dl1 )   ( dl1 )   (  .   .   .   .  A45  .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .   .   .  A54  .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              J = p_KcolA45(ia)
              FF(inode+lofsl1) = FF(inode+lofsl1) &
                              - rvanka%Dmultipliers(4,5)*p_DA45(ia)*p_Dvector(J+ioffsetl2)
              FF(inode+lofsl2) = FF(inode+lofsl2) &
                              - rvanka%Dmultipliers(5,4)*p_DA54(ia)*p_Dvector(J+ioffsetl1)

              ! Whereever we find a DOF that couples to another DOF on the
              ! same element, we put that to both A-blocks of our local matrix.
              do k=1,nnvel
                if (j .eq. IdofGlobal(k)) then
                  AA (inode+lofsl1,k+lofsl2) = p_DA45(ia)*rvanka%Dmultipliers(4,5)
                  AA (inode+lofsl2,k+lofsl1) = p_DA54(ia)*rvanka%Dmultipliers(5,4)
                  exit
                end if
              end do
            end do
          end if

          ! Process A33 if it exists

          if (associated(p_KdiagonalA33)) then

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   I1  .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.

            daux = rvanka%Dmultipliers(3,3)
            FF(1+lofsp) = FF(1+lofsp) &
                        - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
            AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

          end if

          ! Process A66 if it exists

          if (associated(p_KdiagonalA66)) then

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   I2 ) ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.

            daux = rvanka%Dmultipliers(6,6)
            FF(1+lofsxi) = FF(1+lofsxi) &
                        - daux*p_DA66(p_KdiagonalA66(IEL))*p_Dvector(IEL+ioffsetxi)
            AA(1+lofsxi,1+lofsxi) = daux*p_DA66(p_KdiagonalA66(IEL))

          end if

          ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

          ib1=p_KldB(idof)
          ib2=p_KldB(idof+1)-1
          do ib = ib1,ib2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .  B1   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .  B2   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .  B1  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .  B2  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolB(ib)

            ! primal equation
            daux = p_Dvector(j+ioffsetp)
            FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanka%Dmultipliers(1,3)
            FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanka%Dmultipliers(2,3)

            ! dual equation
            daux2 = p_Dvector(j+ioffsetxi)
            FF(inode+lofsl1) = FF(inode+lofsl1)-p_DB1(ib)*daux2 * rvanka%Dmultipliers(4,6)
            FF(inode+lofsl2) = FF(inode+lofsl2)-p_DB2(ib)*daux2 * rvanka%Dmultipliers(5,6)

            ! Do not incorporate the B-matrices into AA yet; this will come later!
          end do

          ! The mass matrix defect.
          if (associated(p_KldM)) then
            ! We assume: multiplier of A(1,4) = multiplier of A(2,5)
            daux = rvanka%Dmultipliers(1,4)

            ! We assume: multiplier of A(4,1) = multiplier of A(5,2)
            daux2 = rvanka%Dmultipliers(4,1)

            ia1 = p_KldM(idof)
            ia2 = p_KldM(idof+1)-1
            do ia = ia1,ia2

              J = p_KcolM(ia)

              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .   .   .  aM   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   (  .   .   .   .  aM   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !   ( dl1 )   ( dl1 )   ( bM   .   .   .   .   .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .  bM   .   .   .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              FF(inode+lofsu) = FF(inode+lofsu)-daux*p_DM(ia)*p_Dvector(J+ioffsetl1)
              FF(inode+lofsv) = FF(inode+lofsv)-daux*p_DM(ia)*p_Dvector(J+ioffsetl2)

              FF(inode+lofsl1) = FF(inode+lofsl1)-daux2*p_DM(ia)*p_Dvector(J+ioffsetu)
              FF(inode+lofsl2) = FF(inode+lofsl2)-daux2*p_DM(ia)*p_Dvector(J+ioffsetv)

              ! Whereever we find a DOF that couples to another DOF on the
              ! same element, we put that to both A-blocks of our local matrix.
              do k=1,nnvel
                if (j .eq. IdofGlobal(k)) then
                  AA (inode+lofsu,k+lofsl1) = daux*p_DM(ia)
                  AA (inode+lofsv,k+lofsl2) = daux*p_DM(ia)

                  AA (inode+lofsl1,k+lofsu) = daux2*p_DM(ia)
                  AA (inode+lofsl2,k+lofsv) = daux2*p_DM(ia)
                  ! AA (k+lofsl1,inode+lofsu) = daux2*p_DM(ia)
                  ! AA (k+lofsl2,inode+lofsv) = daux2*p_DM(ia)
                  exit
                end if
              end do
            end do
          end if

          ! Ok, up to now, all loops are clean and vectoriseable. Now the only
          ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
          ! We have to find in the B-matrices the column that corresponds
          ! to our element and pressure DOF IEL - which makes it necessary
          ! to compare the column numbers in KcolB with IEL.
          ! Remember: The column numbers in B correspond to the pressure-DOF`s
          ! and so to element numbers.
          !
          ! Btw: Each row of B has at most two entries:
          !
          !      IEL                              IEL
          !   |--------|             |--------|--------|
          !   |        |             |        |        |
          !   |   P1   |      or     |   P2   X   P1   |
          !   |        |             |        |        |
          ! --|---X----|--           |--------|--------|
          !
          ! Either two (if the velocity DOF is an edge with two neighbouring
          ! elements) or one (if the velocity DOF is at an edge on the boundary
          ! and there is no neighbour).
          do ib = ib1,ib2

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   ( D1  D2   .   .   .   .  ) ( p  )
            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .  D1  D2   .  ) ( xi )
            !
            ! In AA, we simultaneously set up (locally):
            !
            !   (  .   .  B1   .   .   .  )
            !   (  .   .  B2   .   .   .  )
            !   ( D1  D2   .   .   .   .  )
            !   (  .   .   .   .   .  B1  )
            !   (  .   .   .   .   .  B2  )
            !   (  .   .   .  D1  D2   .  )

            if (p_KcolB(ib) .eq. IEL) then

              J = p_KcolB(ib)

              ! Get the entries in the B-matrices.
              ! Primal equation
              AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanka%Dmultipliers(1,3)
              AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanka%Dmultipliers(2,3)

              ! The same way, get DD1 and DD2.
              ! Note that DDi has exacty the same matrix structrure as BBi and is noted
              ! as 'transposed matrix' only because of the transposed-flag.
              ! So we can use "ib" as index here to access the entry of DDi:
              AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanka%Dmultipliers(3,1)
              AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanka%Dmultipliers(3,2)

              ! The same for the dual equation
              AA(inode+lofsl1,1+lofsxi) = p_DB1(ib) * rvanka%Dmultipliers(4,6)
              AA(inode+lofsl2,1+lofsxi) = p_DB2(ib) * rvanka%Dmultipliers(5,6)

              AA(1+lofsxi,inode+lofsl1) = p_DD1(ib) * rvanka%Dmultipliers(6,4)
              AA(1+lofsxi,inode+lofsl2) = p_DD2(ib) * rvanka%Dmultipliers(6,5)

              ! Build the pressure entry in the local defect vector:
              !   f_i = (f_i-Aui) - D_i pi
              ! or more precisely (as D is roughly B^T):
              !   f_i = (f_i-Aui) - (B^T)_i pi
              FF(1+lofsp) = FF(1+lofsp) &
                          - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                          - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)

              ! The same for the dual pressure
              FF(1+lofsxi) = FF(1+lofsxi) &
                          - AA(1+lofsxi,inode+lofsl1)*p_Dvector(idof+ioffsetl1) &
                          - AA(1+lofsxi,inode+lofsl2)*p_Dvector(idof+ioffsetl2)

              ! Quit the loop - the other possible entry belongs to another
              ! element, not to the current one
              exit
            end if
          end do ! ib

        end do ! inode

        ! Now we make a defect-correction approach for this system:
        !
        !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
        !                                     -----------
        !                                        =d~
        !
        ! Here the 'projection' operator simply converts the small
        ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
        ! of the same size as x - what is easy using the number of
        ! the DOF`s on the element.
        !
        ! For C, we use our local AA, i.e. applying C^{-1} means to
        ! solve the local system AA dd = FF for dd. The local defect dd is then
        ! added back to the global solution vector.

        call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

        if (ilapackInfo .eq. 0) then

          ! Ok, we got the update vector in FF. Incorporate this now into our
          ! solution vector with the update formula
          !
          !  x_{n+1} = x_n + domega * y

          do inode=1,nnvel
            ! Update of the primal velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetu) &
              = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
            p_Dvector(idofGlobal(inode)+ioffsetv) &
              = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)

            ! Update of the dual velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetl1) &
              = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * FF(inode+lofsl1)
            p_Dvector(idofGlobal(inode)+ioffsetl2) &
              = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * FF(inode+lofsl2)
          end do

          p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                    domega * FF(1+lofsp)

          p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + &
                                    domega * FF(1+lofsxi)

        else if (ilapackInfo .lt. 0) then

          call output_line (&
              'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNSSOCQ1TQ0fullCoupCfFB')

        end if

        ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
        ! coarse grid with all boundaries = Dirichlet.
        ! In this case, nothing must be changed in the vector!

      end do ! iel

    ! -------------------------------------------------------------------------
    case (1)
      ! VANKA applied to 9x9 primal system

      do ielidx=1,size(IelementList)

        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)

        ! Clear the 'primal local system matrix'.
        do ib=1,nnprimal
          do ia=1,nnprimal
            AA(ia,ib) = 0.0_DP
          end do
        end do

        ! We now have the element
        !
        ! +---------+                       +----3----+
        ! |         |                       |         |
        ! |   IEL   |   with DOF`s          4    P    2
        ! |         |                       |    Q0   |
        ! +---------+                       +----1----+
        !
        !
        ! Fetch the pressure P on the current element into FF.
        ! The numbers of the DOF`s coincide with the definition
        ! in dofmapping.f90!

        ! Get the primal pressure
        FF(1+lofsp) = p_Drhs(iel+ioffsetp)

        ! Get the velocity DOF`s on the current element.
        ! We assume: DOF 1..4 = edge.
        ! That is the same implementation as in dofmapping.f90!
        IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

        ! Loop over all U-nodes of that element.
        do inode=1,nnvel

          ! Get the DOF we have to tackle:
          idof = IdofGlobal(inode)

          ! Set FF initially to the value of the right hand
          ! side vector that belongs to our current DOF corresponding
          ! to inode.

          ! Primal equation
          FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
          FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

          ! What do we have at this point?
          ! FF     : "local" RHS vector belonging to the DOF`s on the
          !          current element
          ! AA     : Diagonal entries of A belonging to these DOF`s
          !
          ! And at the moment:
          ! idof      : number of current DOF on element IEL
          ! inode     : "local" number of DOF on element IEL, i.e.
          !              number of the edge
          !

          ! Extract those entries in the A-, B- and M-matrices to our local
          ! matrix AA, which belong to DOF`s in our current solution vector.
          !
          ! At first build: fi = fi-Aui

          ia1 = p_KldA11(idof)
          ia2 = p_KldA11(idof+1)-1
          do ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - ( A11  .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .  A22  .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )

            J = p_KcolA11(ia)

            ! Primal equation:
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanka%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanka%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanka%Dmultipliers(1,1)
                AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanka%Dmultipliers(2,2)
                exit
              end if
            end do
          end do

          ! Process the 'off-diagonal' matrices A12 and A21

          if (associated(p_KldA12)) then
            ia1 = p_KldA12(idof)
            ia2 = p_KldA12(idof+1)-1
            do ia = ia1,ia2
              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .  A12  .   .   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   ( A21  .   .   .   .   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !                                                   ( l1 )
              !                                                   ( l2 )
              !                                                   ( xi )

              J = p_KcolA12(ia)
              FF(inode+lofsu) = FF(inode+lofsu) &
                              - rvanka%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
              FF(inode+lofsv) = FF(inode+lofsv) &
                              - rvanka%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)

              ! Whereever we find a DOF that couples to another DOF on the
              ! same element, we put that to both A-blocks of our local matrix.
              do k=1,nnvel
                if (j .eq. IdofGlobal(k)) then
                  AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanka%Dmultipliers(1,2)
                  AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanka%Dmultipliers(2,1)
                  exit
                end if
              end do
            end do
          end if

          ! Process A33 if it exists

          if (associated(p_KdiagonalA33)) then

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   I1  .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.

            daux = rvanka%Dmultipliers(3,3)
            FF(1+lofsp) = FF(1+lofsp) &
                        - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
            AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

          end if

          ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

          ib1=p_KldB(idof)
          ib2=p_KldB(idof+1)-1
          do ib = ib1,ib2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .  B1   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .  B2   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )

            J = p_KcolB(ib)

            ! primal equation
            daux = p_Dvector(j+ioffsetp)
            FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanka%Dmultipliers(1,3)
            FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanka%Dmultipliers(2,3)

            ! Do not incorporate the B-matrices into AA yet; this will come later!
          end do

          ! The mass matrix defect.
          if (associated(p_KldM)) then
            ! We assume: multiplier of A(1,4) = multiplier of A(2,5)
            daux = rvanka%Dmultipliers(1,4)

            ia1 = p_KldM(idof)
            ia2 = p_KldM(idof+1)-1
            do ia = ia1,ia2

              J = p_KcolM(ia)

              ! Calculate:
              !
              !   ( du  ) = ( du  ) - (  .   .   .  aM   .   .  ) ( u  )
              !   ( dv  )   ( dv  )   (  .   .   .   .  aM   .  ) ( v  )
              !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
              !                                                   ( l1 )
              !                                                   ( l2 )
              !                                                   ( xi )

              FF(inode+lofsu) = FF(inode+lofsu)-daux*p_DM(ia)*p_Dvector(J+ioffsetl1)
              FF(inode+lofsv) = FF(inode+lofsv)-daux*p_DM(ia)*p_Dvector(J+ioffsetl2)

              ! Whereever we find a DOF that couples to another DOF on the
              ! same element, we put that to both A-blocks of our local matrix.
              do k=1,nnvel
                if (j .eq. IdofGlobal(k)) then
                  AA (inode+lofsu,k+lofsl1) = daux*p_DM(ia)
                  AA (inode+lofsv,k+lofsl2) = daux*p_DM(ia)
                  exit
                end if
              end do
            end do
          end if

          ! Ok, up to now, all loops are clean and vectoriseable. Now the only
          ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
          ! We have to find in the B-matrices the column that corresponds
          ! to our element and pressure DOF IEL - which makes it necessary
          ! to compare the column numbers in KcolB with IEL.
          ! Remember: The column numbers in B correspond to the pressure-DOF`s
          ! and so to element numbers.
          !
          ! Btw: Each row of B has at most two entries:
          !
          !      IEL                              IEL
          !   |--------|             |--------|--------|
          !   |        |             |        |        |
          !   |   P1   |      or     |   P2   X   P1   |
          !   |        |             |        |        |
          ! --|---X----|--           |--------|--------|
          !
          ! Either two (if the velocity DOF is an edge with two neighbouring
          ! elements) or one (if the velocity DOF is at an edge on the boundary
          ! and there is no neighbour).
          do ib = ib1,ib2

            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
            !   ( dp  )   ( dp  )   ( D1  D2   .   .   .   .  ) ( p  )
            !                                                   ( l1 )
            !                                                   ( l2 )
            !                                                   ( xi )
            !
            ! In AA, we simultaneously set up (locally):
            !
            !   (  .   .  B1   .   .   .  )
            !   (  .   .  B2   .   .   .  )
            !   ( D1  D2   .   .   .   .  )
            !   (  .   .   .   .   .   .  )
            !   (  .   .   .   .   .   .  )
            !   (  .   .   .   .   .   .  )

            if (p_KcolB(ib) .eq. IEL) then

              J = p_KcolB(ib)

              ! Get the entries in the B-matrices.
              ! Primal equation
              AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanka%Dmultipliers(1,3)
              AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanka%Dmultipliers(2,3)

              ! The same way, get DD1 and DD2.
              ! Note that DDi has exacty the same matrix structrure as BBi and is noted
              ! as 'transposed matrix' only because of the transposed-flag.
              ! So we can use "ib" as index here to access the entry of DDi:
              AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanka%Dmultipliers(3,1)
              AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanka%Dmultipliers(3,2)

              ! Build the pressure entry in the local defect vector:
              !   f_i = (f_i-Aui) - D_i pi
              ! or more precisely (as D is roughly B^T):
              !   f_i = (f_i-Aui) - (B^T)_i pi
              FF(1+lofsp) = FF(1+lofsp) &
                          - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                          - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)

              ! Quit the loop - the other possible entry belongs to another
              ! element, not to the current one
              exit
            end if
          end do ! ib

        end do ! inode

        ! Now we make a defect-correction approach for this system:
        !
        !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
        !                                     -----------
        !                                        =d~
        !
        ! Here the 'projection' operator simply converts the small
        ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
        ! of the same size as x - what is easy using the number of
        ! the DOF`s on the element.
        !
        ! For C, we use our local AA, i.e. applying C^{-1} means to
        ! solve the local system AA dd = FF for dd. The local defect dd is then
        ! added back to the global solution vector.
        !
        ! We only set up the primal equation, so we only have to solve a
        ! 9x9 subsystem of the complete 18x18 system. Lapack can handle this
        ! using the 'leading dimension'...

        call DGESV (nnprimal, 1, AA, nnld, Ipiv, FF(1:nnprimal), nnprimal, ilapackInfo)

        if (ilapackInfo .eq. 0) then

          ! Ok, we got the update vector in FF. Incorporate this now into our
          ! solution vector with the update formula
          !
          !  x_{n+1} = x_n + domega * y

          do inode=1,nnvel
            ! Update of the primal velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetu) &
              = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
            p_Dvector(idofGlobal(inode)+ioffsetv) &
              = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)

          end do

          p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                    domega * FF(1+lofsp)

        else if (ilapackInfo .lt. 0) then

          call output_line (&
              'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNSSOCQ1TQ0fullCoupCfFB')

        end if

        ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
        ! coarse grid with all boundaries = Dirichlet.
        ! In this case, nothing must be changed in the vector!

      end do ! iel

    ! -------------------------------------------------------------------------
    case (2)
      ! VANKA applied to 9x9 dual system

      do ielidx=1,size(IelementList)

        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)

        ! Clear the 'dual local system matrix'.
        do ib=nnprimal+1,2*nnprimal
          do ia=nnprimal+1,2*nnprimal
            AA(ia,ib) = 0.0_DP
          end do
        end do

        ! We now have the element
        !
        ! +---------+                       +----3----+
        ! |         |                       |         |
        ! |   IEL   |   with DOF`s          4    P    2
        ! |         |                       |    Q0   |
        ! +---------+                       +----1----+
        !
        !
        ! Fetch the pressure P on the current element into FF.
        ! The numbers of the DOF`s coincide with the definition
        ! in dofmapping.f90!

        ! Get the dual pressure
        FF(1+lofsxi) = p_Drhs(iel+ioffsetxi)

        ! Get the velocity DOF`s on the current element.
        ! We assume: DOF 1..4 = edge.
        ! That is the same implementation as in dofmapping.f90!
        IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)

        ! Loop over all U-nodes of that element.
        do inode=1,nnvel

          ! Get the DOF we have to tackle:
          idof = IdofGlobal(inode)

          ! Set FF initially to the value of the right hand
          ! side vector that belongs to our current DOF corresponding
          ! to inode.

          ! dual equation
          FF(inode+lofsl1) = p_Drhs(idof+ioffsetl1)
          FF(inode+lofsl2) = p_Drhs(idof+ioffsetl2)

          ! What do we have at this point?
          ! FF     : "local" RHS vector belonging to the DOF`s on the
          !          current element
          ! AA     : Diagonal entries of A belonging to these DOF`s
          !
          ! And at the moment:
          ! idof      : number of current DOF on element IEL
          ! inode     : "local" number of DOF on element IEL, i.e.
          !              number of the edge
          !
          ! At first build: fi = fi-Aui

          ia1 = p_KldA11(idof)
          ia2 = p_KldA11(idof+1)-1
          do ia = ia1,ia2
            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .  A44  .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .  A55  .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolA11(ia)

            ! dual equation
            FF(inode+lofsl1) = FF(inode+lofsl1) &
                            - rvanka%Dmultipliers(4,4)*p_DA44(ia)*p_Dvector(J+ioffsetl1)
            FF(inode+lofsl2) = FF(inode+lofsl2) &
                            - rvanka%Dmultipliers(5,5)*p_DA55(ia)*p_Dvector(J+ioffsetl2)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsl1,k+lofsl1) = p_DA44(ia)*rvanka%Dmultipliers(4,4)
                AA (inode+lofsl2,k+lofsl2) = p_DA55(ia)*rvanka%Dmultipliers(5,5)
                exit
              end if
            end do
          end do


          ! Process the 'off-diagonal' matrices A45 and A54 if they exist
          if (associated(p_KldA45)) then
            ia1 = p_KldA45(idof)
            ia2 = p_KldA45(idof+1)-1
            do ia = ia1,ia2
              ! Calculate:
              !
              !                                                   ( u  )
              !                                                   ( v  )
              !                                                   ( p  )
              !   ( dl1 ) = ( dl1 )   (  .   .   .   .  A45  .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .   .   .  A54  .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              J = p_KcolA45(ia)
              FF(inode+lofsl1) = FF(inode+lofsl1) &
                              - rvanka%Dmultipliers(4,5)*p_DA45(ia)*p_Dvector(J+ioffsetl2)
              FF(inode+lofsl2) = FF(inode+lofsl2) &
                              - rvanka%Dmultipliers(5,4)*p_DA54(ia)*p_Dvector(J+ioffsetl1)

              ! Whereever we find a DOF that couples to another DOF on the
              ! same element, we put that to both A-blocks of our local matrix.
              do k=1,nnvel
                if (j .eq. IdofGlobal(k)) then
                  AA (inode+lofsl1,k+lofsl2) = p_DA45(ia)*rvanka%Dmultipliers(4,5)
                  AA (inode+lofsl2,k+lofsl1) = p_DA54(ia)*rvanka%Dmultipliers(5,4)
                  exit
                end if
              end do
            end do
          end if

          ! Process A66 if it exists

          if (associated(p_KdiagonalA66)) then

            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   I2 ) ( xi )
            !
            ! IEL is the pressure DOF which we have to tackle.

            daux = rvanka%Dmultipliers(6,6)
            FF(1+lofsxi) = FF(1+lofsxi) &
                        - daux*p_DA66(p_KdiagonalA66(IEL))*p_Dvector(IEL+ioffsetxi)
            AA(1+lofsxi,1+lofsxi) = daux*p_DA66(p_KdiagonalA66(IEL))

          end if

          ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

          ib1=p_KldB(idof)
          ib2=p_KldB(idof+1)-1
          do ib = ib1,ib2
            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .   .   .  B1  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .  B2  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

            J = p_KcolB(ib)

            ! dual equation
            daux2 = p_Dvector(j+ioffsetxi)
            FF(inode+lofsl1) = FF(inode+lofsl1)-p_DB1(ib)*daux2 * rvanka%Dmultipliers(4,6)
            FF(inode+lofsl2) = FF(inode+lofsl2)-p_DB2(ib)*daux2 * rvanka%Dmultipliers(5,6)

            ! Do not incorporate the B-matrices into AA yet; this will come later!
          end do

          ! The mass matrix defect.
          if (associated(p_KldM)) then
            ! We assume: multiplier of A(4,1) = multiplier of A(5,2)
            daux2 = rvanka%Dmultipliers(4,1)

            ia1 = p_KldM(idof)
            ia2 = p_KldM(idof+1)-1
            do ia = ia1,ia2

              J = p_KcolM(ia)

              ! Calculate:
              !
              !                                                   ( u  )
              !                                                   ( v  )
              !                                                   ( p  )
              !   ( dl1 ) = ( dl1 )   ( bM   .   .   .   .   .  ) ( l1 )
              !   ( dl2 )   ( dl2 )   (  .  bM   .   .   .   .  ) ( l2 )
              !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )

              FF(inode+lofsl1) = FF(inode+lofsl1)-daux2*p_DM(ia)*p_Dvector(J+ioffsetu)
              FF(inode+lofsl2) = FF(inode+lofsl2)-daux2*p_DM(ia)*p_Dvector(J+ioffsetv)

              ! Whereever we find a DOF that couples to another DOF on the
              ! same element, we put that to both A-blocks of our local matrix.
              do k=1,nnvel
                if (j .eq. IdofGlobal(k)) then
                  AA (inode+lofsl1,k+lofsu) = daux2*p_DM(ia)
                  AA (inode+lofsl2,k+lofsv) = daux2*p_DM(ia)
                  exit
                end if
              end do
            end do
          end if

          ! Ok, up to now, all loops are clean and vectoriseable. Now the only
          ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
          ! We have to find in the B-matrices the column that corresponds
          ! to our element and pressure DOF IEL - which makes it necessary
          ! to compare the column numbers in KcolB with IEL.
          ! Remember: The column numbers in B correspond to the pressure-DOF`s
          ! and so to element numbers.
          !
          ! Btw: Each row of B has at most two entries:
          !
          !      IEL                              IEL
          !   |--------|             |--------|--------|
          !   |        |             |        |        |
          !   |   P1   |      or     |   P2   X   P1   |
          !   |        |             |        |        |
          ! --|---X----|--           |--------|--------|
          !
          ! Either two (if the velocity DOF is an edge with two neighbouring
          ! elements) or one (if the velocity DOF is at an edge on the boundary
          ! and there is no neighbour).
          do ib = ib1,ib2

            ! Calculate:
            !
            !                                                   ( u  )
            !                                                   ( v  )
            !                                                   ( p  )
            !   ( dl1 ) = ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
            !   ( dxi )   ( dxi )   (  .   .   .  D1  D2   .  ) ( xi )
            !
            ! In AA, we simultaneously set up (locally):
            !
            !   (  .   .   .   .   .   .  )
            !   (  .   .   .   .   .   .  )
            !   (  .   .   .   .   .   .  )
            !   (  .   .   .   .   .  B1  )
            !   (  .   .   .   .   .  B2  )
            !   (  .   .   .  D1  D2   .  )

            if (p_KcolB(ib) .eq. IEL) then

              J = p_KcolB(ib)

              ! Get the entries in the B-matrices.
              ! Primal equation
              AA(inode+lofsl1,1+lofsxi) = p_DB1(ib) * rvanka%Dmultipliers(4,6)
              AA(inode+lofsl2,1+lofsxi) = p_DB2(ib) * rvanka%Dmultipliers(5,6)

              ! The same way, get DD1 and DD2.
              ! Note that DDi has exacty the same matrix structrure as BBi and is noted
              ! as 'transposed matrix' only because of the transposed-flag.
              ! So we can use "ib" as index here to access the entry of DDi:
              AA(1+lofsxi,inode+lofsl1) = p_DD1(ib) * rvanka%Dmultipliers(6,4)
              AA(1+lofsxi,inode+lofsl2) = p_DD2(ib) * rvanka%Dmultipliers(6,5)

              ! Build the pressure entry in the local defect vector:
              !   f_i = (f_i-Aui) - D_i pi
              ! or more precisely (as D is roughly B^T):
              !   f_i = (f_i-Aui) - (B^T)_i pi
              FF(1+lofsxi) = FF(1+lofsxi) &
                          - AA(1+lofsxi,inode+lofsl1)*p_Dvector(idof+ioffsetl1) &
                          - AA(1+lofsxi,inode+lofsl2)*p_Dvector(idof+ioffsetl2)

              ! Quit the loop - the other possible entry belongs to another
              ! element, not to the current one
              exit
            end if
          end do ! ib

        end do ! inode

        ! Now we make a defect-correction approach for this system:
        !
        !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
        !                                     -----------
        !                                        =d~
        !
        ! Here the 'projection' operator simply converts the small
        ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
        ! of the same size as x - what is easy using the number of
        ! the DOF`s on the element.
        !
        ! For C, we use our local AA, i.e. applying C^{-1} means to
        ! solve the local system AA dd = FF for dd. The local defect dd is then
        ! added back to the global solution vector.
        !
        ! We only set up the primal equation, so we only have to solve a
        ! 9x9 subsystem of the complete 18x18 system. Lapack can handle this
        ! using the 'leading dimension'...
        !
        ! Note that we explicitely using the fact that AA(10+i*18,j) = AA(10,i+j),
        ! so LAPACK will really work on the submatrix A(10:18,10:18) !

        call DGESV (nnprimal, 1, AA(nnprimal+1,nnprimal+1), nnld, &
            Ipiv, FF(nnprimal+1), nnprimal, ilapackInfo)

        if (ilapackInfo .eq. 0) then

          ! Ok, we got the update vector in FF. Incorporate this now into our
          ! solution vector with the update formula
          !
          !  x_{n+1} = x_n + domega * y

          do inode=1,nnvel
            ! Update of the dual velocity vectors
            p_Dvector(idofGlobal(inode)+ioffsetl1) &
              = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * FF(inode+lofsl1)
            p_Dvector(idofGlobal(inode)+ioffsetl2) &
              = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * FF(inode+lofsl2)
          end do

          p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + &
                                    domega * FF(1+lofsxi)

        else if (ilapackInfo .lt. 0) then

          call output_line (&
              'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNSSOCQ1TQ0fullCoupCfFB')

        end if

        ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
        ! coarse grid with all boundaries = Dirichlet.
        ! In this case, nothing must be changed in the vector!

      end do ! iel

    end select

  end subroutine

! *****************************************************************************
! Problem class: VANKA variants for 3D Navier-Stokes problems
! *****************************************************************************

!<subroutine>

  subroutine vanka_init3DNavierStokes (rmatrix,rvanka)

!<description>
  ! Initialises the VANKA variant for 3D Navier-Stokes problems
  ! for conformal discretisations. Checks if the "3D-Navier-Stokes" VANKA variant
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! The substructure rvanka%rvanka3DNavSt is intitialised according
  ! to the information provided in rmatrix.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in), target :: rmatrix
  !</input>

!<inputoutput>
  ! t_vankaPointer3DSPNavSt structure that saves algorithm-specific parameters.
  type(t_vanka), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

    integer :: i,j
    logical :: bextended
    type(t_blockDiscretisation), pointer :: p_rblockDiscr

    bextended = .false.

    ! Matrix must be 4x4.
    if ((rmatrix%nblocksPerCol .ne. 4) .or. (rmatrix%nblocksPerRow .ne. 4)) then
      call output_line ('System matrix is not 4x4.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
      call sys_halt()
    end if

    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(4,:) must be (virtually) transposed. All matrices must be double precision.
    do i=1,4
      do j=1,4

        if (lsysbl_isSubmatrixPresent(rmatrix,i,j)) then

          if (i .le. 3) then
            if (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .ne. 0) then
              call output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
              call sys_halt()
            end if
          else
            if ((i .le. 3) .and. &
               (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .eq. 0)) then
              call output_line ('B1/B2/B3 submatrices must be virtually',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
              call output_line ('transposed (LSYSSC_MSPEC_TRANSPOSED)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
              call sys_halt()
            end if
          end if

          if ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX9)) then
            call output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
            call sys_halt()
          end if

          if (rmatrix%RmatrixBlock(i,j)%cdataType .ne. ST_DOUBLE) then
            call output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
            call sys_halt()
          end if

          ! For scaled matrices, we have to use an extended sub-version of VANKA.
          if ((rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 1.0_DP) .and. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 0.0_DP)) then
            bextended = .true.
          end if

          if ((i .eq. j) .and. &
              (rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 1.0_DP) ) then
            bextended = .true.
          end if

        end if ! neq != 0
      end do
    end do

    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,4)%NA .ne. rmatrix%RmatrixBlock(4,1)%NA) .or. &
        (rmatrix%RmatrixBlock(1,4)%NEQ .ne. rmatrix%RmatrixBlock(4,1)%NCOLS)) then
      call output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,4)%NA .ne. rmatrix%RmatrixBlock(4,2)%NA) .or. &
        (rmatrix%RmatrixBlock(2,4)%NEQ .ne. rmatrix%RmatrixBlock(4,2)%NCOLS)) then
      call output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(3,4)%NA .ne. rmatrix%RmatrixBlock(4,3)%NA) .or. &
        (rmatrix%RmatrixBlock(3,4)%NEQ .ne. rmatrix%RmatrixBlock(4,3)%NCOLS)) then
      call output_line ('Structure of B3 and B3^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
      call sys_halt()
    end if

    ! Fill the output structure with data of the matrices.
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanka%rvanka3DNavSt%p_DA )
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,4),&
        rvanka%rvanka3DNavSt%p_DB1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,4),&
        rvanka%rvanka3DNavSt%p_DB2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,4),&
        rvanka%rvanka3DNavSt%p_DB3)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,1),&
        rvanka%rvanka3DNavSt%p_DD1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,2),&
        rvanka%rvanka3DNavSt%p_DD2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,3),&
        rvanka%rvanka3DNavSt%p_DD3)
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,4),&
        rvanka%rvanka3DNavSt%p_KcolB)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,4), &
        rvanka%rvanka3DNavSt%p_KldB )
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanka%rvanka3DNavSt%p_KcolA)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanka%rvanka3DNavSt%p_KldA )
    if (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanka%rvanka3DNavSt%p_KdiagonalA)
    else
      rvanka%rvanka3DNavSt%p_KdiagonalA => rvanka%rvanka3DNavSt%p_KldA
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,4,4)) then

      ! The matrix must be of format 7 or 9.
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,4),&
          rvanka%rvanka3DNavSt%p_DA44 )

      if (rmatrix%RmatrixBlock(4,4)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(4,4), &
                                rvanka%rvanka3DNavSt%p_KdiagonalA44)
      else
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,4), &
                                rvanka%rvanka3DNavSt%p_KdiagonalA44)
      end if

      ! The presence of A(4,4) forces the extended VANKA to be used
      bextended = .true.

    end if

    if (bextended) rvanka%csubsubtype = 1

    ! What is with A22? Is it the same as A11?
    if (.not. lsyssc_isMatrixContentShared (&
        rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2)) ) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
          rvanka%rvanka3DNavSt%p_DA22 )
    end if

    ! What is with A33? Is it the same as A11?
    if (.not. lsyssc_isMatrixContentShared (&
        rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(3,3)) ) then
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
          rvanka%rvanka3DNavSt%p_DA33 )
    end if

    ! What is with A12 and A21? Do they exist? With a scale factor = 1.0?
!    IF (lsysbl_isSubmatrixPresent(rmatrix,1,2) .AND. &
!        (rmatrix%RmatrixBlock(1,2)%dscaleFactor .EQ. 1.0_DP)) THEN
!      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
!          rvanka%rvanka2DNavSt%p_DA12 )
!
!      CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,2),&
!          rvanka%rvanka2DNavSt%p_KcolA12)
!      CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,2), &
!          rvanka%rvanka2DNavSt%p_KldA12 )
!
!      ! Get the structure. It is assumed that A12 and A21 have the same!
!      IF (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
!        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,2), &
!                                rvanka%rvanka2DNavSt%p_KdiagonalA12)
!      ELSE
!        rvanka%rvanka2DNavSt%p_KdiagonalA12 => rvanka%rvanka2DNavSt%p_KldA12
!      END IF
!
!      IF (.NOT. lsysbl_isSubmatrixPresent(rmatrix,2,1)) THEN
!        CALL output_line ('If A12 is given, A21 must also be given!',&
!            OU_CLASS_ERROR,OU_MODE_STD,'vanka_init2DNavierStokesO')
!        CALL sys_halt()
!      END IF
!
!      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
!          rvanka%rvanka2DNavSt%p_DA21 )
!    END IF

    ! Get the multiplication factors of the submatrices.
    rvanka%rvanka3DNavSt%Dmultipliers(1:4,1:4) = &
        rmatrix%RmatrixBlock(1:4,1:4)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest

    if (.not. associated(p_rblockDiscr)) then
      call output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
      call sys_halt()
    end if

    ! Get the discretisation structure of U,V,W and P from the block
    ! discretisation structure.
    rvanka%rvanka3DNavSt%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanka%rvanka3DNavSt%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanka%rvanka3DNavSt%p_rspatialDiscrW => p_rblockDiscr%RspatialDiscr(3)
    rvanka%rvanka3DNavSt%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(4)

    if ((rvanka%rvanka3DNavSt%p_rspatialDiscrU%inumFESpaces .ne. &
         rvanka%rvanka3DNavSt%p_rspatialDiscrV%inumFESpaces) .or. &
        (rvanka%rvanka3DNavSt%p_rspatialDiscrU%inumFESpaces .ne. &
         rvanka%rvanka3DNavSt%p_rspatialDiscrW%inumFESpaces)) then
      call output_line (&
          'Discretisation structures of X-, Y- and Z-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
      call sys_halt()
    end if

    if ((rvanka%rvanka3DNavSt%p_rspatialDiscrP%inumFESpaces .ne. 1) .and. &
        (rvanka%rvanka3DNavSt%p_rspatialDiscrP%inumFESpaces .ne. &
         rvanka%rvanka3DNavSt%p_rspatialDiscrU%inumFESpaces)) then
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      call output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DNavierStokes')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_3DNavierStokes (rvanka3DNavSt, rvector, rrhs, domega, &
      csubtype, csubsubtype)

!<description>
  ! This routine applies the VANKA variant for 3D Navier-Stokes problems
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANKA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector

  ! The subtype of VANKA that should handle the above problem class.
  ! One of the VANKATP_xxxx constants, e.g. VANKATP_DIAGONAL.
  integer :: csubtype

  ! The sub-subtype of VANKA that should handle the above problem class.
  ! =0: use standard VANKA. =1: use extended VANKA (e.g. with different
  !     multipliers in the matrices)
  integer :: csubsubtype

!</input>

!<inputoutput>
  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointer3DNavSt), intent(inout) :: rvanka3DNavSt
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ielementdist
    integer, dimension(:), pointer :: p_IelementList
    type(t_elementDistribution), pointer :: p_relementDistrU
    type(t_elementDistribution), pointer :: p_relementDistrV
    type(t_elementDistribution), pointer :: p_relementDistrW
    type(t_elementDistribution), pointer :: p_relementDistrP

    ! 3D Navier Stokes problem.

    ! Loop through the element distributions of the velocity.
    do ielementdist = 1,rvanka3DNavSt%p_rspatialDiscrU%inumFESpaces

      ! Get the corresponding element distributions of U, V, W and P.
      p_relementDistrU => &
          rvanka3DNavSt%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanka3DNavSt%p_rspatialDiscrV%RelementDistr(ielementdist)
      p_relementDistrW => &
          rvanka3DNavSt%p_rspatialDiscrW%RelementDistr(ielementdist)

      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      if (rvanka3DNavSt%p_rspatialDiscrP%inumFESpaces .gt. 1) then
        p_relementDistrP => &
            rvanka3DNavSt%p_rspatialDiscrP%RelementDistr(ielementdist)
      else
        p_relementDistrP => &
            rvanka3DNavSt%p_rspatialDiscrP%RelementDistr(1)
      end if

      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y- and Z-velocity (and to that
      ! of the pressure).
      call storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)

      ! Which element combination do we have now?
      if ((elem_getPrimaryElement(p_relementDistrU%celement) .eq. EL_Q1T_3D) .and. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .eq. EL_Q1T_3D) .and. &
          (elem_getPrimaryElement(p_relementDistrW%celement) .eq. EL_Q1T_3D) .and. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .eq. EL_Q0_3D)) then
        ! Q1~/Q1~/Q1~/Q0 discretisation

        ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
        select case (csubtype)
        case (VANKATP_DIAGONAL)
!          ! Diagonal VANKA; check if A12 exists.
!          IF (.NOT. ASSOCIATED(rvanka3DNavSt%p_DA12)) THEN
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            if (rvanka3DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! Uniform discretisation
              call vanka_3DSPQ1TQ0simple (rvanka3DNavSt, &
                  rvector, rrhs, domega)
            else
              ! Conformal discretisation
              call vanka_3DSPQ1TQ0simpleConf (rvanka3DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
!          ELSE
!            ! Apply the conformal VANKA that allows different matrices
!            ! in A11, A12, A21 and A22!
!            CALL vanka_2DSPQ1TQ0simpleCoupConf (rvanka3DNavSt, &
!                  rvector, rrhs, domega,p_IelementList)
!          END IF

        case (VANKATP_FULL)
!          ! Full VANKA; check if A12 exists.
!          IF (.NOT. ASSOCIATED(rvanka3DNavSt%p_DA12)) THEN
            ! Call the VANKA subsolver to apply VANKA to our current element list.
            ! Note: Up to now, there is no 'full' variant -- has to be implemented!
            if (rvanka3DNavSt%p_rspatialDiscrU%inumFESpaces .eq. 1) then
              ! uniform discretisation;
              ! here, use the same as for the general conformal discretisation.
              ! Could be speeded up by introducing another variant...
              call vanka_3DSPQ1TQ0fullConf (rvanka3DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            else
              ! Conformal discretisation
              call vanka_3DSPQ1TQ0fullConf (rvanka3DNavSt, &
                  rvector, rrhs, domega,p_IelementList)
            end if
!          ELSE
!            ! Apply the conformal VANKA that allows different matrices
!            ! in A11, A12, A21 and A22!
!            ! If we have multiplication factors, we even have to use an extended
!            ! version of this.
!            IF (csubsubtype .EQ. 0) THEN
!              CALL vanka_2DSPQ1TQ0fullCoupConf (rvanka3DNavSt, &
!                  rvector, rrhs, domega,p_IelementList)
!            ELSE
!              CALL vanka_2DNSQ1TQ0fullCoupConfExt (rvanka3DNavSt, &
!                  rvector, rrhs, domega,p_IelementList)
!            END IF
!          END IF

        case default
          call output_line (&
              'Unknown VANKA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_3DNavierStokes')
          call sys_halt()

        end select

      else
        call output_line (&
            'Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_3DNavierStokes')
        call sys_halt()

      end if

    end do

  end subroutine

  ! ***************************************************************************
  ! 3D Navier-Stokes VANKA, simple diagonal version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    / A              B1 \
  !    |      A         B2 |
  !    |           A    B3 |
  !    \ D1^T D2^T D3^T 0  /
  !
  ! with D1/D2/D3 having the same structure as B1/B2/D3 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1, B2 and B3 are the same matrices as D1, D2 and D3. The only
  ! difference: Some rows in B1/B2/B3 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_init3DSPQ1TQ0simple (rmatrix,rvanka)

!<description>
  ! Checks if the "3D-Navier-Stokes-Q1T-Q0" VANKA variant can be applied to
  ! the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>

!<output>
  ! t_vankaPointer3DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer3DNavSt), intent(out) :: rvanka
!</output>

!</subroutine>

    integer :: i,j

    ! Matrix must be 4x4.
    if ((rmatrix%nblocksPerCol .ne. 4) .or. (rmatrix%nblocksPerRow .ne. 4)) then
      call output_line ('System matrix is not 4x4.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
      call sys_halt()
    end if

    ! A(1:2,1:3) must not be virtually transposed and of format 9.
    ! A(4,:) must be (virtually) transposed. All matrices must be double precision.
    do i=1,4
      do j=1,4

        if (lsysbl_isSubmatrixPresent(rmatrix,i,j)) then

          if (i .le. 3) then
            if (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .ne. 0) then
              call output_line ('Transposed submatrices not supported.',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
              call sys_halt()
            end if
          else
            if (iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) &
                .eq. 0) then
              call output_line ('B1/B2/B3 submatrices must be virtually '//&
                  'transposed (LSYSSC_MSPEC_TRANSPOSED)',&
                  OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
              call sys_halt()
            end if
          end if

          if ((rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
              (rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX9)) then
            call output_line ('Only format 7 and 9 matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
            call sys_halt()
          end if

          if (rmatrix%RmatrixBlock(i,j)%cdataType .ne. ST_DOUBLE) then
            call output_line ('Only double precision matrices supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
            call sys_halt()
          end if

          if (rmatrix%RmatrixBlock(i,j)%dscaleFactor .ne. 1.0_DP) then
            call output_line ('Scaled matrices not supported.',&
                OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
            call sys_halt()
          end if

        end if ! neq != 0
      end do
    end do

    ! The structure of A(1,4) must be identical to A(4,1) and
    ! that of A(2,4) must be identical to A(4,2).
    if ((rmatrix%RmatrixBlock(1,4)%NA .ne. rmatrix%RmatrixBlock(4,1)%NA) .or. &
        (rmatrix%RmatrixBlock(1,4)%NEQ .ne. rmatrix%RmatrixBlock(4,1)%NCOLS)) then
      call output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,4)%NA .ne. rmatrix%RmatrixBlock(4,2)%NA) .or. &
        (rmatrix%RmatrixBlock(2,4)%NEQ .ne. rmatrix%RmatrixBlock(4,2)%NCOLS)) then
      call output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(3,4)%NA .ne. rmatrix%RmatrixBlock(4,3)%NA) .or. &
        (rmatrix%RmatrixBlock(3,4)%NEQ .ne. rmatrix%RmatrixBlock(4,3)%NCOLS)) then
      call output_line ('Structure of B3 and B3^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init3DSPQ1TQ0simple')
      call sys_halt()
    end if

    ! Fill the output structure with data of the matrices.
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),rvanka%p_DA )
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,4),rvanka%p_DB1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,4),rvanka%p_DB2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,4),rvanka%p_DB3)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,1),rvanka%p_DD1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,2),rvanka%p_DD2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,3),rvanka%p_DD3)
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,4),rvanka%p_KcolB)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,4), rvanka%p_KldB )
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),rvanka%p_KcolA)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), rvanka%p_KldA )
    if (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                               rvanka%p_KdiagonalA)
    else
      rvanka%p_KdiagonalA => rvanka%p_KldA
    end if

    ! Get the multiplication factors of the submatrices
    rvanka%Dmultipliers(:,:) = rmatrix%RmatrixBlock(1:4,1:4)%dscaleFactor

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_3DSPQ1TQ0simple (rvanka, rvector, rrhs, domega)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 3D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointer3DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer3DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DB3
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2
    real(DP), dimension(:), pointer :: p_DD3

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IfacesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetw,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 6
    integer, parameter :: lofsw = 12
    integer, parameter :: lofsp = 18
    real(DP) :: daux
    !REAL(DP), DIMENSION(6,6) :: Dmult

    ! Local arrays for informations about one element
    real(DP), dimension(6) :: AA,BB1,BB2,BB3,DD1,DD2,DD3
    real(DP), dimension(19) :: FF,UU
    integer, dimension(6) :: idofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.
    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DB3 => rvanka%p_DB3
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2
    p_DD3 => rvanka%p_DD3

    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanka%Dmultipliers(:,:)

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd, 3rd and 4th component of the solution
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetw = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetp = ioffsetw+rvector%RvectorBlock(3)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================
    !
    ! If you need a description with fancy ASCII arts abou the algorithm
    ! that is implemented below, scroll up to the 2D version.

    ! So we start with a loop over all elements
    do iel=1,NEL

      ! Fetch the pressure P on the current element into FFP
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 6 U-nodes of that element.
      do inode=1,6

        ! Set idof to the DOF that belongs to our face node:
        idof = p_IfacesAtElement(inode,iel)

        ! Write the number of the face node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))

        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANKA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let us hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        FF(inode+lofsw) = p_Drhs(idof+ioffsetw)

        ! At first build: fi = fi-Aui
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)

          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          !FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
          FF(inode+lofsw) = FF(inode+lofsw)-p_DB3(ib)*daux

          ! For support of scaled matrices, use the following lines; currently switched off.
          !FF(inode)       = FF(inode)      -Dmult(1,4)*p_DB1(ib)*daux
          !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,4)*p_DB2(ib)*daux
          !FF(inode+lofsw) = FF(inode+lofsw)-Dmult(3,4)*p_DB3(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local BX and DX.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        ! Either two (if the velocity DOF is an face with two neighbouring
        ! elements) or one (if the velocity DOF is at an face on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            BB3(inode) = p_DB3(ib)

            ! The same way, get DD1, DD2 and DD3.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            DD3(inode) = p_DD3(ib)

            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,4)*p_DB1(ib)
            !BB2(inode) = Dmult(2,4)*p_DB2(ib)
            !BB3(inode) = Dmult(3,4)*p_DB3(ib)
            !DD1(inode) = Dmult(4,1)*p_DD1(ib)
            !DD2(inode) = Dmult(4,2)*p_DD2(ib)
            !DD3(inode) = Dmult(4,3)*p_DD3(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv) &
                        - DD3(inode)*p_Dvector(idof+ioffsetw)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_3DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,BB3,DD1,DD2,DD3)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,6
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
        p_Dvector(idofGlobal(inode)+ioffsetw) &
          = p_Dvector(idofGlobal(inode)+ioffsetw) + domega * UU(inode+lofsw)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_3DSPQ1TQ0simpleConf (rvanka, rvector, rrhs, domega,IelementList)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 3D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_3DSPQ1TQ0simpleConf is the same as vanka_3DSPQ1TQ0simple except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer3DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer3DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DB3
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2
    real(DP), dimension(:), pointer :: p_DD3

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IfacesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetw,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 6
    integer, parameter :: lofsw = 12
    integer, parameter :: lofsp = 18
    real(DP) :: daux

    ! Local arrays for informations about one element
    real(DP), dimension(6) :: AA,BB1,BB2,BB3,DD1,DD2,DD3
    real(DP), dimension(19) :: FF,UU
    integer, dimension(6) :: idofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.
    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DB3 => rvanka%p_DB3
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2
    p_DD3 => rvanka%p_DD3

    ! For support of scaled matrices, use the following line; currently switched off.
    !Dmult(:,:) = rvanka%Dmultipliers(:,:)

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd, 3rd and 4th component of the solution
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetw = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetp = ioffsetw+rvector%RvectorBlock(3)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================
    !
    ! If you need a description with fancy ASCII arts abou the algorithm
    ! that is implemented below, scroll up to the 2D version.

    ! So we start with a loop over all elements
    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Fetch the pressure P on the current element into FFP
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 6 U-nodes of that element.
      do inode=1,6

        ! Set idof to the DOF that belongs to our face node:
        idof = p_IfacesAtElement(inode,iel)

        ! Write the number of the face node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))

        ! For support of scaled matrices, use the following line; currently switched off.
        ! Node that this way, VANKA would not support a different scaling factor for
        ! A(1,1) than for A(2,2)! Let us hope that this is nowhere used!
        !AA(inode) = Dmult(1,1)*p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        FF(inode+lofsw) = p_Drhs(idof+ioffsetw)

        ! At first build: fi = fi-Aui
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)

          ! For support of scaled matrices, use the following line; currently switched off.
          !daux = Dmult(1,1)*p_DA(ia)
          !FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          !FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          !FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
          FF(inode+lofsw) = FF(inode+lofsw)-p_DB3(ib)*daux

          ! For support of scaled matrices, use the following lines; currently switched off.
          !FF(inode)       = FF(inode)      -Dmult(1,4)*p_DB1(ib)*daux
          !FF(inode+lofsv) = FF(inode+lofsv)-Dmult(2,4)*p_DB2(ib)*daux
          !FF(inode+lofsw) = FF(inode+lofsw)-Dmult(3,4)*p_DB3(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local BX and DX.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        ! Either two (if the velocity DOF is an face with two neighbouring
        ! elements) or one (if the velocity DOF is at an face on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)
            BB3(inode) = p_DB3(ib)

            ! The same way, get DD1, DD2 and DD3.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)
            DD3(inode) = p_DD3(ib)

            ! For support of scaled matrices, use the following lines; currently switched off.
            !BB1(inode) = Dmult(1,4)*p_DB1(ib)
            !BB2(inode) = Dmult(2,4)*p_DB2(ib)
            !BB3(inode) = Dmult(3,4)*p_DB3(ib)
            !DD1(inode) = Dmult(4,1)*p_DD1(ib)
            !DD2(inode) = Dmult(4,2)*p_DD2(ib)
            !DD3(inode) = Dmult(4,3)*p_DD3(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv) &
                        - DD3(inode)*p_Dvector(idof+ioffsetw)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_3DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,BB3,DD1,DD2,DD3)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,6
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
        p_Dvector(idofGlobal(inode)+ioffsetw) &
          = p_Dvector(idofGlobal(inode)+ioffsetw) + domega * UU(inode+lofsw)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)

    end do ! iel

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine vanka_getcorr_3DSPQ1TQ0simple (UU,FF,AA,BB1,BB2,BB3,DD1,DD2,DD3)

!<description>
  ! This routine solves a 19x19 Jacobi-type Schur complement system for three
  ! velocity vectors and one pressure vector. It is used as auxiliary
  ! routine in the simple VANKA solver to calculate an update vector
  ! for velocity/pressure.
!</description>

!<input>
  ! Diagonal elements of the local system matrix.
  real(DP), dimension(6), intent(in) :: AA

  ! Entries in the submatrix B1.
  real(DP), dimension(6), intent(in) :: BB1

  ! Entries in the submatrix B2.
  real(DP), dimension(6), intent(in) :: BB2

  ! Entries in the submatrix B3.
  real(DP), dimension(6), intent(in) :: BB3

  ! Entries in the submatrix D1.
  real(DP), dimension(6), intent(in) :: DD1

  ! Entries in the submatrix D2.
  real(DP), dimension(6), intent(in) :: DD2

  ! Entries in the submatrix D3.
  real(DP), dimension(6), intent(in) :: DD3

  ! Local RHS vector; FF(1..6)=X-velocity, FF(7..12)=Y-velocity,
  ! FF(13..18)=Y-velocity, FF(19)=pressure.
  real(DP), dimension(19), intent(in) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..6)=X-velocity, UU(7..12)=Y-velocity,
  ! UU(13..18)=Y-velocity, UU(19)=pressure.
  real(DP), dimension(19), intent(out) :: UU
!</output>

!</subroutine>

    ! local variables

    integer :: inode
    real(DP) :: PP,dpres
    real(DP), dimension(19) :: AI,dff

    integer, parameter :: lofsv = 6
    integer, parameter :: lofsw = 12
    integer, parameter :: lofsp = 18

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial!
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    do inode=1,6

      ! Quick check if everything is ok - we do not want to divide by 0.
      if (AA(inode)*AA(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      ! AI(.) saves the diagonal matrix S^-1:
      AI(inode)=1.0_DP/AA(inode)

    end do

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 6x1 matrix
    ! b) S^-1 is a diagonal matrix, given by the 6 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = 0.0_DP
    dff = FF

    do inode = 1,6
      dpres        = dpres - AI(inode)*( &
                     DD1(inode)*BB1(inode)+ &
                     DD2(inode)*BB2(inode)+ &
                     DD3(inode)*BB3(inode))
      dff(1+lofsp) = dff(1+lofsp) - AI(inode)*( &
                     DD1(inode)*dff(inode)+ &
                     DD2(inode)*dff(inode+lofsv)+ &
                     DD3(inode)*dff(inode+lofsw))
    end do

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF`s are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    if (dpres*dpres .lt. 1E-20_DP)  then
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      return
    endif

    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:

    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP

    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:

    do inode=1,6
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
      UU(inode+lofsw) = AI(inode)*(dff(inode+lofsw)-BB3(inode)*PP)
    end do

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine vanka_getcorr_3DSPQ1TQ0simple2 (UU,FF,AA1,AA2,AA3,&
                                          BB1,BB2,BB3,DD1,DD2,DD3,di)

!<description>
  ! This routine solves a 19x19 Jacobi-type Schur complement system for three
  ! velocity vectors and one pressure vector. It is used as auxiliary
  ! routine in the simple VANKA solver to calculate an update vector
  ! for velocity/pressure.
  !
  ! In contrast to vanka_getcorr_3DSPQ1TQ0simple, the three diagonal blocks
  ! in the matrix may be different from each other.
!</description>

!<input>
  ! Diagonal elements of the local system matrix A11
  real(DP), dimension(6), intent(in) :: AA1

  ! Diagonal elements of the local system matrix A22
  real(DP), dimension(6), intent(in) :: AA2

  ! Diagonal elements of the local system matrix A33
  real(DP), dimension(6), intent(in) :: AA3

  ! Entries in the submatrix B1.
  real(DP), dimension(6), intent(in) :: BB1

  ! Entries in the submatrix B2.
  real(DP), dimension(6), intent(in) :: BB2

  ! Entries in the submatrix B3.
  real(DP), dimension(6), intent(in) :: BB3

  ! Entries in the submatrix D1.
  real(DP), dimension(6), intent(in) :: DD1

  ! Entries in the submatrix D2.
  real(DP), dimension(6), intent(in) :: DD2

  ! Entries in the submatrix D3.
  real(DP), dimension(6), intent(in) :: DD3

  ! Entry in the submatrix I (usually =0).
  ! This is the matrix in the diagonal block of the pressure, which is usually
  ! zero in saddle point problems.
  real(DP), intent(in) :: di

  ! Local RHS vector; FF(1..6)=X-velocity, FF(7..12)=Y-velocity,
  ! FF(13..18)=Y-velocity, FF(19)=pressure.
  real(DP), dimension(19), intent(in) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..6)=X-velocity, UU(7..12)=Y-velocity,
  ! UU(13..18)=Y-velocity, UU(19)=pressure.
  real(DP), dimension(19), intent(out) :: UU
!</output>

!</subroutine>

    ! local variables

    integer :: inode
    real(DP) :: PP,dpres
    real(DP), dimension(19) :: AI1,AI2,AI3,dff

    integer, parameter :: lofsv = 6
    integer, parameter :: lofsw = 12
    integer, parameter :: lofsp = 18

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T i )
    !      ( D1  D2  I  )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t I ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !       D^t u + Ip = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !      D^t S^-1 (f-Bp) + Ip = g
    !
    !      <=>  Ip -D^t S^-1 B p  =  g - D^t S^-1 f
    !              ***********       **************
    !                 =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial!
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    do inode=1,6

      ! Quick check if everything is ok - we do not want to divide by 0.
      if (AA1(inode)*AA1(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      if (AA2(inode)*AA2(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      if (AA3(inode)*AA3(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      ! AI(.) saves the diagonal matrix S^-1:
      AI1(inode)=1.0_DP/AA1(inode)
      AI2(inode)=1.0_DP/AA2(inode)
      AI3(inode)=1.0_DP/AA3(inode)

    end do

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !      ( I - D^t S^-1 B ) p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 6x1 matrix
    ! b) S^-1 is a diagonal matrix, given by the 6 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           = (I - D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = di
    dff = FF

    do inode = 1,6
      dpres        = dpres &
                   - AI1(inode)*DD1(inode)*BB1(inode) &
                   - AI2(inode)*DD2(inode)*BB2(inode) &
                   - AI3(inode)*DD3(inode)*BB3(inode)

      dff(1+lofsp) = dff(1+lofsp) &
                   - AI1(inode)*DD1(inode)*dff(inode) &
                   - AI2(inode)*DD2(inode)*dff(inode+lofsv) &
                   - AI3(inode)*DD3(inode)*dff(inode+lofsw)
    end do

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF`s are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    if (dpres*dpres .lt. 1E-20_DP)  then
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      return
    endif

    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:

    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP

    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:

    do inode=1,6
      UU(inode)       = AI1(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI2(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
      UU(inode+lofsw) = AI3(inode)*(dff(inode+lofsw)-BB3(inode)*PP)
    end do

  end subroutine

  ! ***************************************************************************
  ! 3D Navier-Stokes VANKA, 'full' version.
  ! Supports Q1~/Q0 discretisation only.
  ! Matrix must be of the form
  !
  !    / A              B1 \
  !    |      A         B2 |
  !    |           A    B3 |
  !    \ D1^T D2^T D3^T 0  /
  !
  ! with D1/D2/D3 having the same structure as B1/B2/B3 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1, B2 and B3 are the same matrices as D1, D2 and D3. The only
  ! difference: Some rows in B1/B2/B3 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

!<subroutine>

  subroutine vanka_3DSPQ1TQ0fullConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 3D Navier-Stokes problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_3DSPQ2QP1fullConf is the same as vanka_3DSPQ2QP1full except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer3DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer3DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DB3
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2
    real(DP), dimension(:), pointer :: p_DD3

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IfacesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 6      ! Q1T = 6 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! QQ0 = 1 DOF`s per pressure
    integer, parameter :: nnld = 3*nnvel+nnpressure
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetw,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsw = 2*nnvel
    integer, parameter :: lofsp = 3*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.
    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DB3 => rvanka%p_DB3
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2
    p_DD3 => rvanka%p_DD3

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetw = ioffsetv+rvector%RvectorBlock(2)%NEQ
    ioffsetp = ioffsetw+rvector%RvectorBlock(3)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! So we start with a loop over all elements in the list
    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..6 = face.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:6) = p_IfacesAtElement(1:6,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode
        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
        FF(inode+lofsw) = p_Drhs(idof+ioffsetw)

        ! At first build: fi = fi-Aui
        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
          FF(inode+lofsw) = FF(inode+lofsw)-daux*p_Dvector(J+ioffsetw)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to all A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = daux
              AA (inode+lofsv,k+lofsv) = daux
              AA (inode+lofsw,k+lofsw) = daux
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
          FF(inode+lofsw) = FF(inode+lofsw)-p_DB3(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        ! Either two (if the velocity DOF is an face with two neighbouring
        ! elements) or one (if the velocity DOF is at an face on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            AA(inode,      lofsp+1) = p_DB1(ib)
            AA(inode+lofsv,lofsp+1) = p_DB2(ib)
            AA(inode+lofsw,lofsp+1) = p_DB3(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(lofsp+1,inode)       = p_DD1(ib)
            AA(lofsp+1,inode+lofsv) = p_DD2(ib)
            AA(lofsp+1,inode+lofsw) = p_DD3(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(lofsp+1,inode)*p_Dvector(idof) &
                        - AA(lofsp+1,inode+lofsv)*p_Dvector(idof+ioffsetv) &
                        - AA(lofsp+1,inode+lofsw)*p_Dvector(idof+ioffsetw)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  .............. :: :: :: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(4,4) :: :: :: )
      !     (                            AA( 5, 5) .............. :: :: :: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA( 8, 8) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
          p_Dvector(idofGlobal(inode)+ioffsetw) &
            = p_Dvector(idofGlobal(inode)+ioffsetw) + domega * FF(inode+lofsw)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_3DSPQ1TQ0fullConf')

      end if
      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************
  ! 2D Navier-Stokes VANKA, simple diagonal and full version.
  ! Supports P1~/P0 discretisation.
  !
  ! with D1/D2 having the same structure as B1/B2 and the 'transposed'
  ! flag set (LSYSSC_MSPEC_TRANSPOSED).
  ! In general, B1 and B2 are the same matrices as D1 and D2. The only
  ! difference: Some rows in B1/B2 may be replaced by zero lines to implement
  ! Dirichlet boundary conditions.
  ! ***************************************************************************

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPP1TP0simpleConf (rvanka, rvector, rrhs, domega,IelementList)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with P1~/P0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ1TQ0simpleConf is the same as vanka_2DSPQ1TQ0simple except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 3
    integer, parameter :: lofsp = 6
    real(DP) :: daux

    ! Local arrays for informations about one element
    real(DP), dimension(3) :: AA,BB1,BB2,DD1,DD2
    real(DP), dimension(7) :: FF,UU
    integer, dimension(3) :: IdofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      do inode=1,3

        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the four velocity unknowns and the
        ! pressure unknown on the current element; these five unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity
        ! components and B being an (2*4) x 1 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |        |             |        |        |
        !   --|---X----|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                     BB1(1) )
      !     (        AA(2)                              BB1(2) )
      !     (               AA(3)                       BB1(3) )
      !     (                      AA(1)                BB2(1) )
      !     (                             AA(2)         BB2(2) )
      !     (                                    AA(3)  BB2(3) )
      !     ( DD1(1) DD1(2) DD1(3) DD2(1) DD2(2) DD2(3)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_2DSPP1TP0simple2 (UU,FF,AA,AA,BB1,BB2,DD1,DD2,0.0_DP)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,3
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)

    end do ! iel

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine vanka_getcorr_2DSPP1TP0simple2 (UU,FF,AA1,AA2,BB1,BB2,DD1,DD2,di)

!<description>
  ! This routine solves a 9x9 Jacobi-type Schur complement system for two
  ! velocity vectors and one pressure vector. It is used as auxiliary
  ! routine in the simple VANKA solver to calculate an update vector
  ! for velocity/pressure.
  !
  ! In contrast to vanka_getcorr_2DSPQ1TQ0simple, the two diagonal blocks
  ! in the matrix may be different from each other.
!</description>

!<input>
  ! Diagonal elements of the local system matrix A11
  real(DP), dimension(3), intent(in) :: AA1

  ! Diagonal elements of the local system matrix A22
  real(DP), dimension(3), intent(in) :: AA2

  ! Entries in the submatrix B1.
  real(DP), dimension(3), intent(in) :: BB1

  ! Entries in the submatrix B2.
  real(DP), dimension(3), intent(in) :: BB2

  ! Entries in the submatrix D1.
  real(DP), dimension(3), intent(in) :: DD1

  ! Entries in the submatrix D2.
  real(DP), dimension(3), intent(in) :: DD2

  ! Entry in the submatrix I (usually =0).
  ! This is the matrix in the diagonal block of the pressure, which is usually
  ! zero in saddle point problems.
  real(DP), intent(in) :: di

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  real(DP), dimension(7), intent(in) :: FF
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  real(DP), dimension(7), intent(out) :: UU
!</output>

!</subroutine>

    ! local variables

    integer :: inode
    real(DP) :: PP,dpres
    real(DP), dimension(3) :: AI1,AI2
    real(DP), dimension(7) :: dff

    integer, parameter :: lofsv = 3
    integer, parameter :: lofsp = 6

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA1(1)                                    BB1(1) )
    !     (        AA1(2)                             BB1(2) )
    !     (               AA1(3)                      BB1(3) )
    !     (                      AA2(1)               BB2(1) )
    !     (                             AA2(2)        BB2(2) )
    !     (                                    AA2(3) BB2(3) )
    !     ( DD1(1) DD1(2) DD1(3) DD2(1) DD2(2) DD2(3) di1    )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T i )
    !      ( D1  D2  I  )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t I ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !       D^t u + Ip = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !      D^t S^-1 (f-Bp) + Ip = g
    !
    !      <=>  Ip -D^t S^-1 B p  =  g - D^t S^-1 f
    !              ***********       **************
    !                 =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial!
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    do inode=1,3

      ! Quick check if everything is ok - we do not want to divide by 0.
      if (AA1(inode)*AA1(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      if (AA2(inode)*AA2(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      ! AI(.) saves the diagonal matrix S^-1:

      AI1(inode)=1E0_DP/AA1(inode)
      AI2(inode)=1E0_DP/AA2(inode)

    end do

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !      ( I - D^t S^-1 B ) p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 3x1 matrix
    ! b) S^-1 is a diagonal matrix, given by the 3 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           = (I - D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = di
    dff = FF

    do inode = 1,3
      dpres        = dpres &
                   - AI1(inode)*DD1(inode)*BB1(inode) &
                   - AI2(inode)*DD2(inode)*BB2(inode)
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI1(inode)*DD1(inode)*dff(inode) &
                   - AI2(inode)*DD2(inode)*dff(inode+lofsv)
    end do

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF`s are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    if (dpres*dpres .lt. 1E-20_DP)  then
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      return
    endif

    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:

    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP

    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:

    do inode=1,3
      UU(inode)       = AI1(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI2(inode)*(dff(inode+lofsv)-BB2(inode)*PP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPP1TP0simpleCoupConf (rvanka, rvector, rrhs, domega,IelementList)

!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with P1~/P0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ1TQ0simpleCoupConf is the same as vanka_2DSPQ1TQ0simpleConf,
  ! but supports fully coupled velocity submatrices.
  ! The matrices A11 and A22 must have the same structure. The matrices A12
  ! and A21 must also have the same structure. The structure of A11 and A12
  ! may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12,p_KdiagonalA12
    real(DP), dimension(:), pointer :: p_DA,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp
    integer :: ia1,ia2,ib1,ib2,ia,ib,j
    integer, parameter :: lofsv = 3
    integer, parameter :: lofsp = 6
    real(DP) :: daux

    ! Local arrays for informations about one element
    real(DP), dimension(2*lofsv) :: AA,BB1,BB2,DD1,DD2
    real(DP), dimension(7) :: FF,UU
    integer, dimension(lofsv) :: idofGlobal

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_DA22 => rvanka%p_DA22

    ! Structure of A12 is assumed to be the same as A21
    p_KcolA12 => rvanka%p_KcolA12
    p_KldA12 => rvanka%p_KldA12
    p_KdiagonalA12 => rvanka%p_KdiagonalA12
    p_DA12 => rvanka%p_DA12
    p_DA21 => rvanka%p_DA21

    ! Structure of B1 is assumed to be the same as B2
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        |---X---|                |---X---|
    !        |       |                |       |
    !    old X   1   X old   -->  new X   1   X new
    !        |       |                |       |
    !        |---X---|                |---X---|
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new     new
    !        |---X---|---X---|             |---X---|---X---|
    !        |       |       |             |       |       |
    !    new X   1   X   2   X old --> new X   1   X   2   X new
    !        |       |new    |             |       |newer  |
    !        |---X---|---X---|             |---X---|---X---|
    !           new     old                   new     new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! We now have the element
      !                                      U3/V3
      ! |---------|                       |----X----|
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s    U4/V4 X    P    X U2/V2
      ! |         |                       |         |
      ! |---------|                       |----X----|
      !                                      U1/V1
      !
      ! Fetch the pressure P on the current element into FFP

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Loop over all 4 U-nodes of that element.
      do inode=1,lofsv

        ! Set idof to the DOF that belongs to our edge inode:
        idof = p_IedgesAtElement(inode,iel)

        ! Write the number of the edge/node to idofGlobal:
        idofGlobal(inode) = idof

        ! Put on AA(.) the diagonal entry of matrix A
        AA(inode) = p_DA(p_KdiagonalA(idof))
        AA(inode+lofsv) = p_DA22(p_KdiagonalA(idof))

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the four velocity unknowns and the
        ! pressure unknown on the current element; these five unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 4 x (2*NVT) matrix for the two velocity
        ! components and B being an (2*4) x 1 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |        |             |        |        |
        !   --|---X----|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to
        ! two 4 x 1 submatrices (originally, every velocity couples with
        ! two pressure elements on the neighbour cell, so we have
        ! 2 columns in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)
        end do

        ! Tackle 'offdiagonal' matrices A12 and A21

        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            BB1(inode) = p_DB1(ib)
            BB2(inode) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            DD1(inode) = p_DD1(ib)
            DD2(inode) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - DD1(inode)*p_Dvector(idof) &
                        - DD2(inode)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1)                                                   BB1(1) )
      !     (        AA(2)                                            BB1(2) )
      !     (               AA(3)                                     BB1(3) )
      !     (                      AA(4)                              BB1(4) )
      !     (                             AA(1)                       BB2(1) )
      !     (                                    AA(2)                BB2(2) )
      !     (                                           AA(3)         BB2(3) )
      !     (                                                  AA(4)  BB2(4) )
      !     ( DD1(1) DD1(2) DD1(3) DD1(4) DD2(1) DD2(2) DD2(3) DD2(4)        )
      !
      ! We could theoretically pass this to LAPACK or so to invert it - but
      ! as this is a saddle-point system, we can much faster 'solve' the equation
      ! y := C^{-1} d~ by applying a Schur complement approach. For this purpose,
      ! call the element update routine that calculates the update vector y.

      call vanka_getcorr_2DSPP1TP0CPsimple (UU,FF,AA,BB1,BB2,DD1,DD2,0.0_DP)

      ! Ok, we got the update vector UU. Incorporate this now into our
      ! solution vector with the update formula
      !
      !  x_{n+1} = x_n + domega * y!

      do inode=1,lofsv
        p_Dvector(idofGlobal(inode)) &
          = p_Dvector(idofGlobal(inode)) + domega * UU(inode)
        p_Dvector(idofGlobal(inode)+ioffsetv) &
          = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * UU(inode+lofsv)
      end do

      p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + domega * UU(1+lofsp)

    end do ! iel

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine vanka_getcorr_2DSPP1TP0CPsimple (UU,FF,AA,BB1,BB2,DD1,DD2,di)

!<description>
  ! This routine solves a 7x7 Jacobi-type Schur complement system for two
  ! velocity vectors and one pressure vector. It is used as auxiliary
  ! routine in the simple VANKA solver to calculate an update vector
  ! for velocity/pressure for system where the velocity is fully coupled.
!</description>

!<input>
  ! Diagonal elements of the local system matrix. (Blocks A11 and A22)
  real(DP), dimension(2*3), intent(in) :: AA

  ! Entries in the submatrix B1.
  real(DP), dimension(3), intent(in) :: BB1

  ! Entries in the submatrix B2.
  real(DP), dimension(3), intent(in) :: BB2

  ! Entries in the submatrix D1.
  real(DP), dimension(3), intent(in) :: DD1

  ! Entries in the submatrix D2.
  real(DP), dimension(3), intent(in) :: DD2

  ! Local RHS vector; FF(1..4)=X-velocity, FF(5..8)=Y-velocity,
  ! FF(9)=pressure.
  real(DP), dimension(7), intent(in) :: FF

  ! Entry in the submatrix I (usually =0).
  ! This is the matrix in the diagonal block of the pressure, which is usually
  ! zero in saddle point problems.
  real(DP), intent(in) :: di
!</input>

!<output>
  ! Update vector u with Cu=FF. UU(1..4)=X-velocity, UU(5..8)=Y-velocity,
  ! UU(9)=pressure.
  real(DP), dimension(7), intent(out) :: UU
!</output>

!</subroutine>

    ! local variables

    integer :: inode
    real(DP) :: PP,dpres
    real(DP), dimension(7) :: AI,dff

    integer, parameter :: lofsv = 3
    integer, parameter :: lofsp = 6

    ! This routine uses a Schur-complement approach to solve the
    ! system Cu=FF with
    !
    ! C = ( AA(1)                                     BB1(1) )
    !     (        AA(2)                              BB1(2) )
    !     (               AA(3)                       BB1(3) )
    !     (                      AA(4)                BB2(1) )
    !     (                             AA(5)         BB2(2) )
    !     (                                    AA(6)  BB2(3) )
    !     ( DD1(1) DD1(2) DD1(3) DD2(1) DD2(2) DD2(3)        )
    !
    !   =: ( A       B1 )  =:  ( S   B )
    !      (     A   B2 )      ( D^T 0 )
    !      ( D1  D2     )
    !
    ! What we want to calculate are two things: 1.) a new pressure and
    ! 2.) a new velocity. Both can be calculated from the
    ! RHS using the Schur-Complement approach.
    !
    ! Assume we have a system:
    !
    !  [ S   B ] (u) = (f)
    !  [ D^t 0 ] (p)   (g)
    !
    ! We can write:
    !
    !                u = S^-1 (f-Bp)
    !            D^t u = g
    !
    ! Inserting the first equation into the second one gives:
    !
    !           D^t S^-1 (f-Bp) = g
    !
    !      <=>   -D^t S^-1 B p  =  g - D^t S^-1 f
    !            ***********       **************
    !               =: DP              =: FF(pressure)
    !
    ! Note that DP is a 1x1-system, i.e. a scalar! Therefore
    ! calculating DP^-1 to get p=DP^-1*FF(pressure) is trivial!
    ! So FF(pressure)/DP will be the pressure on the element IEL.
    !
    ! Calculating an update for the velocity
    !
    !      u = S^-1 (f-Bp)
    !
    ! is then also trivial as S (and thus S^-1) is a diagonal matrix.
    !
    ! Here it goes...

    do inode=1,2*lofsv

      ! Quick check if everything is ok - we do not want to divide by 0.
      if (AA(inode)*AA(inode) .lt. 1E-20_DP) then
        ! Set the update vector to 0, cancel.
        UU = 0.0_DP
        return
      end if

      ! AI(.) saves the diagonal matrix S^-1:

      AI(inode)=1E0_DP/AA(inode)

    end do

    ! Factorization loop
    !
    ! What we at first want to calculate is p with
    !
    !         - D^t S^-1 B p  =  g - D^t S^-1 f
    !
    ! To calculate that for local B, S, f and p, consider at first
    ! the dimensions in this system:
    !
    ! a) B is a 3x1 matrix
    ! b) S^-1 is a diagonal matrix, given by the 3 diagonal entries of A
    ! c) D^t S^-1 B is therefore a 1x1 matrix, thus a scalar
    !
    ! So in the factorization loop we can calculate:
    !
    !   DP           =   - (D^T S^-1 B)
    !   FF(pressure) = g - (D^T S^-1 f)
    !
    ! As S and S^-1 are a diagonal matrices, we can exploit
    ! B^T S^-1  =  S^-1 B^T  which saves some multiplications...

    dpres = di
    dff = FF

    do inode = 1,lofsv
      dpres        = dpres &
                   - AI(inode)  *(DD1(inode)*BB1(inode)) &
                   - AI(inode+lofsv)*(DD2(inode)*BB2(inode))
      dff(1+lofsp) = dff(1+lofsp) &
                   - AI(inode)  *(DD1(inode)*dff(inode)) &
                   - AI(inode+lofsv)*(DD2(inode)*dff(inode+lofsv))
    end do

    ! Solution "loop"
    !
    ! Check that DP exists. It may be e.g. ~0 if all velocity DOF`s are Dirichlet
    ! nodes, which implies B=0 and thus leads to DP=0!
    ! (Happens inside of fictitious boundary objects e.g.)

    if (dpres*dpres .lt. 1E-20_DP)  then
      ! Set the update vector to 0, cancel.
      UU = 0.0_DP
      return
    endif

    ! At first we calculate the pressure on element IEL,
    ! which is simply given by multiplying FFP with the
    ! inverte "matrix" DP, i.e.:

    PP          = dff(1+lofsp)/dpres
    UU(1+lofsp) = PP

    ! With the help of the pressure, calculate the velocity.
    ! This can be done again by the Schur-Complement approach using
    !
    !       u = S^-1 (f-Bp)
    !
    ! locally on the current cell:

    do inode=1,lofsv
      UU(inode)       = AI(inode)*(dff(inode)-BB1(inode)*PP)
      UU(inode+lofsv) = AI(inode+lofsv)*(dff(inode+lofsv)-BB2(inode)*PP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPP1TP0fullConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier-Stokes problems with P1~/P0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ2QP1fullConf is the same as vanka_2DSPQ2QP1full except
  ! for IelementList. This parameter allows to specify a subset of elements
  ! where VANKA should be applied to. Other elements than in this list are
  ! ignored!
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 3      ! P1T = 3 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! P0 = 1 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          4    P    2
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..3 = edge.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:nnvel) = p_IedgesAtElement(1:nnvel,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 8 x 2*(NMT) matrix for the two velocity
        ! components and B~ being an (2*3) x 2 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! 2*1 pressure elements on the adjacent cells, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with two square
        ! submatrices A~ of size 3 x 3. The 6 x 2-matrix B~ reduces to
        ! two 3 x 2 submatrices (originally, every velocity couples with
        ! the pressure DOF`s on that cell, so we have
        ! 1 column in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          daux = p_DA(ia)
          FF(inode)       = FF(inode)      -daux*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-daux*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = daux
              AA (inode+nnvel,k+nnvel) = daux
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            AA(inode,      2*nnvel+1) = p_DB1(ib)
            AA(inode+nnvel,2*nnvel+1) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(2*nnvel+1,inode)       = p_DD1(ib)
            AA(2*nnvel+1,inode+nnvel) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(2*nnvel+1,inode)*p_Dvector(idof) &
                        - AA(2*nnvel+1,inode+nnvel)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  .............. :: :: :: )
      !     (    :                  :                                   :AA :: )
      !     (    :                  :                                   :(B1): )
      !     (    ................ AA(3,3) :: :: :: )
      !     (                            AA( 4, 4) .............. :: :: :: )
      !     (                                :                  :       :AA :: )
      !     (                                :                  :       :(B2): )
      !     (                                ............... AA( 6, 6) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)
      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPP1TP0fullConf')

      end if
      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DSPP1TP0fullCoupConf (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier-Stokes problems with P1~/P0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DSPQ1TQ0fullCoupConf is the same as vanka_2DSPQ1TQ0fullConf,
  ! but supports fully coupled velocity submatrices.
  ! The matrices A11 and A22 must have the same structure. The matrices A12
  ! and A21 must also have the same structure. The structure of A11 and A12
  ! may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA
    integer, dimension(:), pointer :: p_KldA,p_KdiagonalA
    real(DP), dimension(:), pointer :: p_DA,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12,p_KdiagonalA12
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 3      ! P1T = 3 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! P0 = 1 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel+nnpressure   ! Q2/Q2/P1 = 9+9+3 = 21 DOF`s per element
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! offset information in arrays
    integer :: ioffsetv,ioffsetp,j
    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA => rvanka%p_KcolA
    p_KldA => rvanka%p_KldA
    p_KdiagonalA => rvanka%p_KdiagonalA
    p_DA => rvanka%p_DA
    p_DA22 => rvanka%p_DA22

    ! Structure of A12 is assumed to be the same as A21
    p_KcolA12 => rvanka%p_KcolA12
    p_KldA12 => rvanka%p_KldA12
    p_KdiagonalA12 => rvanka%p_KdiagonalA12
    p_DA12 => rvanka%p_DA12
    p_DA21 => rvanka%p_DA21

    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          4    P    2
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !
      !
      ! Fetch the pressure P on the current element into FFP.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..3 = edge.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:nnvel) = p_IedgesAtElement(1:nnvel,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode

        FF(inode)       = p_Drhs(idof)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u) = (f)
        !    [ B^t 0 ] (p)   (g)
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in a rectangular
        ! system of the form
        !
        !    [ === A^ === B~ ] (|) = (f1)
        !    [ B~^t       0  ] (u)   (f2)
        !                      (|)   (g )
        !                      (p)
        !
        ! with A^ being an 8 x 2*(NMT) matrix for the two velocity
        ! components and B~ being an (2*3) x 2 matrix that couples the
        ! velocities to the pressure on our current element.
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! 2*1 pressure elements on the adjacent cells, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL!
        !
        !   (d1)  = (f1) - [ === A^ === B~ ] (u1)
        !   (d2)    (f2)   [ B~^t       0  ] (u2)
        !   (dp)    (g )                     (p)
        !
        !
        ! That way, A^ is reduced to a square matrix with four square
        ! submatrices A~ of size 3 x 3. The 6 x 2-matrix B~ reduces to
        ! two 3 x 2 submatrices (originally, every velocity couples with
        ! the pressure DOF`s on that cell, so we have
        ! 1 column in the B-matrix).
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA(idof)
        ia2 = p_KldA(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA(ia)
          FF(inode)       = FF(inode)      -p_DA(ia)*p_Dvector(J)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k) = p_DA(ia)
              AA (inode+nnvel,k+nnvel) = p_DA22(ia)
              exit
            end if
          end do
        end do

        ! Process the 'off-diagonal' matrices A12 and A21

        ia1 = p_KldA12(idof)
        ia2 = p_KldA12(idof+1)-1
        do ia = ia1,ia2
          J = p_KcolA12(ia)
          FF(inode)       = FF(inode)      -p_DA12(ia)*p_Dvector(J+ioffsetv)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DA21(ia)*p_Dvector(J)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode,k+nnvel) = p_DA12(ia)
              AA (inode+nnvel,k) = p_DA21(ia)
              exit
            end if
          end do
        end do

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          J = p_KcolB(ib)
          daux = p_Dvector(j+ioffsetp)
          FF(inode)       = FF(inode)      -p_DB1(ib)*daux
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2
          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices
            AA(inode,      2*nnvel+1) = p_DB1(ib)
            AA(inode+nnvel,2*nnvel+1) = p_DB2(ib)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(2*nnvel+1,inode)       = p_DD1(ib)
            AA(2*nnvel+1,inode+nnvel) = p_DD2(ib)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(2*nnvel+1,inode)*p_Dvector(idof) &
                        - AA(2*nnvel+1,inode+nnvel)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! The only question now will be: What is C^{-1}?
      !
      ! Well, here we have different choices.
      ! For full linear systems, one would choose C=A, which ist the
      ! theoretically best preconditioner. A more simple preconditioner
      ! is a kind of Jacobi-preconditioner, which extracts the main diagonal
      ! entries of A and those lines of the B/D-matrices that correspond
      ! to the DOF`s on the current element. We already set up the preconditioner
      ! in the above variables. It has the form:
      !
      ! C = ( AA(1,1)  ..............    AA( 1, 4) .............. :: :: :: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B1): )
      !     (    ................ AA(3,3)    ............... AA( 4, 6) :: :: :: )
      !     ( AA(4,1)  ..............    AA( 4, 4) .............. :: :: :: )
      !     (    :                  :        :                  :       :AA :: )
      !     (    :                  :        :                  :       :(B2): )
      !     (    ................ AA(6,3)    ............... AA( 6, 6) :: :: :: )
      !     ( ===== AA (D1-block) =====  ======= AA (D2-block) ======          )
      !
      ! To solve this (a little bit larger) system, we invoke LAPACK.

      !CALL DGETRF( nnld, nnld, AA, nnld, Ipiv, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRF) LU decomposition'
      !CALL DGETRS('N', nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo )
      !IF(ilapackInfo.ne.0) PRINT *,'ERROR: LAPACK(DGETRS) back substitution'

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y!

        do inode=1,nnvel
          p_Dvector(idofGlobal(inode)) &
            = p_Dvector(idofGlobal(inode)) + domega * FF(inode)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)

      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DSPQ1TQ0fullCoupConf')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_2DNSP1TP0fullCoupConfExt (rvanka, rvector, rrhs, domega, IelementList)

!<description>
  ! This routine applies the specialised full local system VANKA algorithm for
  ! 2D Navier Stokes optimal control problems with Q1~/Q0 discretisation
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! vanka_2DNSSQ1TQ0fullCoupConf supports fully coupled velocity submatrices.
  ! The matrices A11, A22 must have the same structure.
  ! The matrices A12 and A21 must have the same structure.
  ! The structure of A11 and A12 may be different from each other.
!</description>

!<input>
  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointer2DNavSt), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:) :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,ielidx
    integer :: inode,idof

    integer, dimension(:), pointer :: p_KcolA11
    integer, dimension(:), pointer :: p_KldA11
    real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22
    integer, dimension(:), pointer :: p_KcolA12
    integer, dimension(:), pointer :: p_KldA12
    integer, dimension(:), pointer :: p_KcolB
    integer, dimension(:), pointer :: p_KldB
    real(DP), dimension(:), pointer :: p_DB1
    real(DP), dimension(:), pointer :: p_DB2
    real(DP), dimension(:), pointer :: p_DD1
    real(DP), dimension(:), pointer :: p_DD2
    real(DP), dimension(:), pointer :: p_Da33
    integer, dimension(:), pointer :: p_KdiagonalA33

    ! Triangulation information
    integer :: NEL
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector

    ! Local arrays for informations about one element
    integer, parameter :: nnvel = 3      ! P1T = 3 DOF`s per velocity
    integer, parameter :: nnpressure = 1 ! P0 = 1 DOF`s per pressure
    integer, parameter :: nnld = 2*nnvel + nnpressure
    integer, dimension(nnvel) :: IdofGlobal
    real(DP), dimension(nnld,nnld) :: AA
    real(DP), dimension(nnld) :: FF

    ! Offsets of the 'local' solution parts in the 'local' solution vector
    integer, parameter :: lofsu = 0
    integer, parameter :: lofsv = nnvel
    integer, parameter :: lofsp = 2*nnvel

    ! LAPACK temporary space
    integer :: Ipiv(nnld),ilapackInfo

    ! Offset information in arrays.
    integer :: ioffsetu,ioffsetv,ioffsetp,j

    integer :: ia1,ia2,ib1,ib2,ia,ib,k
    real(DP) :: daux

    ! Get pointers to the system matrix, so we do not have to write
    ! so much - and it is probably faster.

    ! Structure of A11 is assumed to be the same as A22
    p_KcolA11 => rvanka%p_KcolA
    p_KldA11 => rvanka%p_KldA
    p_DA11 => rvanka%p_DA
    p_DA22 => rvanka%p_DA22
    if (.not. associated(p_DA22)) p_DA22 => p_DA11

    ! Structure of A12 is assumed to be the same as A21.
    ! Get A12 and A21 -- except for if the multipliers are =0, then
    ! we switch them off by nullifying the pointers.
    if (rvanka%Dmultipliers(1,2) .ne. 0.0_DP) then
      p_KcolA12 => rvanka%p_KcolA12
      p_KldA12 => rvanka%p_KldA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    else
      nullify(p_KcolA12)
      nullify(p_KldA12)
      nullify(p_DA12 )
      nullify(p_DA21 )
    end if

    p_KcolB => rvanka%p_KcolB
    p_KldB => rvanka%p_KldB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    ! Diagonal submatrices A33 and A66 (if they exist)
    if (rvanka%Dmultipliers(3,3) .ne. 0.0_DP) then
      p_Da33 => rvanka%p_DA33
      p_KdiagonalA33 => rvanka%p_KdiagonalA33
    else
      nullify(p_Da33)
      nullify(p_KdiagonalA33)
    end if

    ! Get pointers to the vectors, RHS, get triangulation information
    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lsysbl_getbase_double (rvector,p_Dvector)
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get the relative offsets of the 2nd and 3rd solution of the component
    ioffsetu = 0
    ioffsetv = rvector%RvectorBlock(1)%NEQ
    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ

    !=======================================================================
    !     Block Gauss-Seidel on Schur Complement
    !=======================================================================

    ! Basic algorithm:
    !
    ! What are we doing here? Well, we want to perform
    ! *preconditioning*, i.e. we have to solve the problem
    !
    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
    !
    ! for a "special" preconditioner C which we define in a moment.
    ! This is equivalent to solving the system
    !
    !   C (x_new)  = x_old
    !
    ! C should be some approximation to A. Imagine our global system:
    !
    !     [ A   B ] (u) = (f)
    !     [ B^t 0 ] (p)   (g)
    !
    ! In the Navier-Stokes equations with (u,p) being the preconditioned
    ! vector, there should be g=0 - but this cannot be assumed
    ! as it does not happen in general.
    ! Now the algorithm for generating a new (u,p) vector from the old
    ! one reads roughly as follows:
    !
    ! a) Restrict to a small part of the domain, in our case to one cell.
    ! b) Fetch all the data (velocity, pressure) on that cell. On the
    !    first cell, we have only "old" velocity entries. These values
    !    are updated and then the calculation proceeds with the 2nd cell.
    !
    !           old                      new
    !        +---X---+                +---X---+
    !        |       |                |       |
    !    old X       X       -->  new X   X   X new
    !        |   1   |                |   1   |
    !        +---X---+                +---X---+
    !           old                      new
    !
    !    From the second cell on, there might be "old" data and "new"
    !    data on that cell - the old data that has not been updated and
    !    perhaps some already updated velocity data from a neighbor cell.
    !
    !           new     old                   new       new
    !        +---X---+---X---+             +---X---+---X---+
    !        |     1 |     2 |             |     1 |     2 |
    !    new X       X       X old --> new X       X       X new
    !        |       |new    |             |       |newer  |
    !        +---X---+---X---+             +---X---+---X---+
    !           new     old                   new       new
    !
    !    These values are updated and then the calculation proceeds
    !    with the next cell.
    !    As can be seen in the above picture, the "new" node in the
    !    middle is even going to be a "newer" node when handled again
    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
    !    Information is updated subsequently by using "old" data and
    !    "new" data from a previous calculation.
    !
    ! So we start with a loop over all elements in the list

    do ielidx=1,size(IelementList)

      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)

      ! Clear the 'local system matrix'.
      AA(:,:) = 0.0_DP

      ! We now have the element
      !
      ! +---------+                       +----3----+
      ! |         |                       |         |
      ! |   IEL   |   with DOF`s          4    P    2
      ! |         |                       |    Q0   |
      ! +---------+                       +----1----+
      !
      !
      ! Fetch the pressure P on the current element into FF.
      ! The numbers of the DOF`s coincide with the definition
      ! in dofmapping.f90!

      ! Get the pressure
      FF(1+lofsp) = p_Drhs(iel+ioffsetp)

      ! Get the velocity DOF`s on the current element.
      ! We assume: DOF 1..3 = edge.
      ! That is the same implementation as in dofmapping.f90!
      IdofGlobal(1:nnvel) = p_IedgesAtElement(1:nnvel,iel)

      ! Loop over all U-nodes of that element.
      do inode=1,nnvel

        ! Get the DOF we have to tackle:
        idof = IdofGlobal(inode)

        ! Set FF initially to the value of the right hand
        ! side vector that belongs to our current DOF corresponding
        ! to inode.

        ! Primal equation
        FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)

        ! What do we have at this point?
        ! FF     : "local" RHS vector belonging to the DOF`s on the
        !          current element
        ! AA     : Diagonal entries of A belonging to these DOF`s
        !
        ! And at the moment:
        ! idof      : number of current DOF on element IEL
        ! inode     : "local" number of DOF on element IEL, i.e.
        !              number of the edge
        !
        ! Now comes the crucial point with the "update": How to
        ! subsequently update the vertex values, such that the whole
        ! thing still converges to the solution, even if a node
        ! is updated more than once? Here, we use a typical
        ! matrix-decomposition approach:
        !
        ! Again consider the problem:
        !
        !    [ A   B ] (u)  = (f  )
        !    [ B^t 0 ] (p)    (g  )
        !
        ! We assume, that all components in the vector (u,p) are
        ! given - except for the velocity and pressure unknowns
        ! on the current element; these 21 unknowns
        ! are located anywhere in the (u,p) vector. The idea is to
        ! shift "everything known" to the right hand side to obtain
        ! a system for only these unknowns!
        !
        ! Extracting all the lines of the system that correspond to
        ! DOF`s on our single element IEL results in rectangular
        ! systems of the form
        !
        !    [ === A^ === B~  ] (| ) = (f1 )
        !    [ B~^t       I1~ ] (u )   (f2 )
        !                       (| )   (g  )
        !                       (p )
        !
        !
        ! B~ is a 8 x 2 matrix: As every velocity couples with at most
        ! 2*1 pressure elements on the adjacent cells, so we have
        ! 2 columns in the B-matrix.
        !
        !        IEL                              IEL
        !     |--------|             |--------|--------|
        !     |        |             |        |        |
        !     |   P    |      or     |   Q    X   P    |
        !     |   X    |             |        |        |
        !   --|--------|--           |--------|--------|
        !
        !
        ! Now, throw all summands to the RHS vector to build a local
        ! 'defect' on our single element IEL.
        !
        !  (d1 ) = (f1 ) -  [ === A^ === B~  ] (| )
        !  (d2 )   (f2 )    [ B~^t       I1~ ] (u )
        !  (dg )   (g  )                       (| )
        !                                      (p )
        !
        ! Extract those entries in the A-, B- and M-matrices to our local
        ! matrix AA, which belong to DOF`s in our current solution vector.
        !
        ! At first build: fi = fi-Aui

        ia1 = p_KldA11(idof)
        ia2 = p_KldA11(idof+1)-1
        do ia = ia1,ia2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - ( A11  .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .  A22  .   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   ) ( p  )

          J = p_KcolA11(ia)

          ! Primal equation:
          FF(inode+lofsu) = FF(inode+lofsu) &
                          - rvanka%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
          FF(inode+lofsv) = FF(inode+lofsv) &
                          - rvanka%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)

          ! Whereever we find a DOF that couples to another DOF on the
          ! same element, we put that to both A-blocks of our local matrix.
          do k=1,nnvel
            if (j .eq. IdofGlobal(k)) then
              AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanka%Dmultipliers(1,1)
              AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanka%Dmultipliers(2,2)
              exit
            end if
          end do
        end do

        ! Process the 'off-diagonal' matrices A12 and A21

        if (associated(p_KldA12)) then
          ia1 = p_KldA12(idof)
          ia2 = p_KldA12(idof+1)-1
          do ia = ia1,ia2
            ! Calculate:
            !
            !   ( du  ) = ( du  ) - (  .  A12  .  ) ( u  )
            !   ( dv  )   ( dv  )   ( A21  .   .  ) ( v  )
            !   ( dp  )   ( dp  )   (  .   .   .  ) ( p  )

            J = p_KcolA12(ia)
            FF(inode+lofsu) = FF(inode+lofsu) &
                            - rvanka%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
            FF(inode+lofsv) = FF(inode+lofsv) &
                            - rvanka%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)

            ! Whereever we find a DOF that couples to another DOF on the
            ! same element, we put that to both A-blocks of our local matrix.
            do k=1,nnvel
              if (j .eq. IdofGlobal(k)) then
                AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanka%Dmultipliers(1,2)
                AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanka%Dmultipliers(2,1)
                exit
              end if
            end do
          end do
        end if

        ! Process A33 if it exists

        if (associated(p_KdiagonalA33)) then

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   I1  ) ( p  )
          !
          ! IEL is the pressure DOF which we have to tackle.

          daux = rvanka%Dmultipliers(3,3)
          FF(1+lofsp) = FF(1+lofsp) &
                      - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
          AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))

        end if

        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi

        ib1=p_KldB(idof)
        ib2=p_KldB(idof+1)-1
        do ib = ib1,ib2
          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .  B1   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .  B2   ) ( v  )
          !   ( dp  )   ( dp  )   (  .   .   .   ) ( p  )

          J = p_KcolB(ib)

          daux = p_Dvector(j+ioffsetp)
          FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanka%Dmultipliers(1,3)
          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanka%Dmultipliers(2,3)

          ! Do not incorporate the B-matrices into AA yet; this will come later!
        end do

        ! Ok, up to now, all loops are clean and vectoriseable. Now the only
        ! somehow 'unclean' loop to determine the local B1, B2, D1 and D2.
        ! We have to find in the B-matrices the column that corresponds
        ! to our element and pressure DOF IEL - which makes it necessary
        ! to compare the column numbers in KcolB with IEL.
        ! Remember: The column numbers in B correspond to the pressure-DOF`s
        ! and so to element numbers.
        !
        ! Btw: Each row of B has at most two entries:
        !
        !      IEL                              IEL
        !   |--------|             |--------|--------|
        !   |        |             |        |        |
        !   |   P1   |      or     |   P2   X   P1   |
        !   |        |             |        |        |
        ! --|---X----|--           |--------|--------|
        !
        ! Either two (if the velocity DOF is an edge with two neighbouring
        ! elements) or one (if the velocity DOF is at an edge on the boundary
        ! and there is no neighbour).
        do ib = ib1,ib2

          ! Calculate:
          !
          !   ( du  ) = ( du  ) - (  .   .   .   ) ( u  )
          !   ( dv  )   ( dv  )   (  .   .   .   ) ( v  )
          !   ( dp  )   ( dp  )   ( D1  D2   .   ) ( p  )
          !
          ! In AA, we simultaneously set up (locally):
          !
          !   (  .   .  B1   )
          !   (  .   .  B2   )
          !   ( D1  D2   .   )

          if (p_KcolB(ib) .eq. IEL) then

            J = p_KcolB(ib)

            ! Get the entries in the B-matrices.
            ! Primal equation
            AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanka%Dmultipliers(1,3)
            AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanka%Dmultipliers(2,3)

            ! The same way, get DD1 and DD2.
            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
            ! as 'transposed matrix' only because of the transposed-flag.
            ! So we can use "ib" as index here to access the entry of DDi:
            AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanka%Dmultipliers(3,1)
            AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanka%Dmultipliers(3,2)

            ! Build the pressure entry in the local defect vector:
            !   f_i = (f_i-Aui) - D_i pi
            ! or more precisely (as D is roughly B^T):
            !   f_i = (f_i-Aui) - (B^T)_i pi
            FF(1+lofsp) = FF(1+lofsp) &
                        - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
                        - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)

            ! Quit the loop - the other possible entry belongs to another
            ! element, not to the current one
            exit
          end if
        end do ! ib

      end do ! inode

      ! Now we make a defect-correction approach for this system:
      !
      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
      !                                     -----------
      !                                        =d~
      !
      ! Here the 'projection' operator simply converts the small
      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
      ! of the same size as x - what is easy using the number of
      ! the DOF`s on the element.
      !
      ! For C, we use our local AA, i.e. applying C^{-1} means to
      ! solve the local system AA dd = FF for dd. The local defect dd is then
      ! added back to the global solution vector.

      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)

      if (ilapackInfo .eq. 0) then

        ! Ok, we got the update vector in FF. Incorporate this now into our
        ! solution vector with the update formula
        !
        !  x_{n+1} = x_n + domega * y

        do inode=1,nnvel
          ! Update of the primal velocity vectors
          p_Dvector(idofGlobal(inode)+ioffsetu) &
            = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
          p_Dvector(idofGlobal(inode)+ioffsetv) &
            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
        end do

        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
                                  domega * FF(1+lofsp)

      else if (ilapackInfo .lt. 0) then

        call output_line (&
            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNSSQ1TQ0fullCoupConf')

      end if

      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
      ! coarse grid with all boundaries = Dirichlet.
      ! In this case, nothing must be changed in the vector!

    end do ! iel

  end subroutine

end module
