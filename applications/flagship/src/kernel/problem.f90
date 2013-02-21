!##############################################################################
!# ****************************************************************************
!# <name> problem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures and subroutines for handling
!#  the complete problem configuration both in space and time.
!#
!# This problem structure is implemented in the most general way in
!# order to support a simple Poisson problem in 1D as well as a
!# multiphysics problem which consists of several subproblems, whereby
!# the spatial dimension of subproblems may by different. As an
!# example consider a fluid-structure interaction problem which
!# consists of a fluid phase "Fluid" and a structure phase
!# "Structure". Each phase represents a problem per se usings its own
!# triangulation, boundary parametrisation, etc. From the programming
!# point of view, we can assume that there exist two black-box
!# solution algorithms for fluid and structure phase,
!# respectively. However, the coupling of both phases cannot be done
!# in a black-box fashion so that the user needs to take care of the
!# coupling in the application code. This type of subproblems will be
!# referred to as "standalone subproblems".
!#
!# As an example consider the following vertical problem hierarchy
!#
!# Example 1:          [FSI]
!#                    /     \
!#             [Fluid]       [Structure]
!#
!# Moreover let the fluid phase be modeled by the compressible Euler
!# equations for two gas phases, that is, subproblem "Fluid" itself is
!# a multi-phase problem itself. Let us further assume that this
!# multi-phase problem uses a common triangulation, boundary
!# parametrisation, etc. That is, each phase is modeled by the
!# single-phase compressible Euler equations and, say, the two phases
!# are described by using the level-set approach. Thus, another scalar
!# transport equation is required to model the fluid phase.
!#
!# Example 2:          [Fluid]  =  [compEuler1]+[compEuler2]+[Transport]
!#
!# In contrast to the above example, we want to treat the fluid
!# problem in a fully coupled way. That is, we want to generate a
!# single system matrix, a single solution vector and a single
!# right-hand side vector. This is possible since the three
!# subproblems are based on the same triangulation, boundary
!# parametrisation, etc. Even further, all subproblems make use of the
!# same hierarchy of meshes (see below). Thus, was refer to
!# subproblems of this kind as "embedded subproblems". Example 2 may
!# also be regarded as flat of horizontal problem hierarchy.
!#
!# However, the aim is to reuse as much as possible from the
!# application code written for, say, the single-phase compressible
!# Euler equations and the scalar transport model. Thus, we must be
!# able to identify subproblem [Transport] which is embedded into
!# problem [Fluid]. Thus, subproblem [Transport] has its "own" problem
!# structure but its content is shared with the encompassing problem
!# structure of the [Fluid] problem.
!# 
!# Technically, the two different types of subproblems are realised as
!# follows. Each (sub)problem is represented by a t_problem structure
!# which has a unique name, say, "FSI".
!# The "FSI" problem is just a meta-problem which means that it does
!# not provide any real data. In contrast, it serves as virtual
!# container for the two standalone subproblems "Fluid" and
!# "Structure". The pointers p_rproblemFirst and p_rproblemLast refer
!# to the first and last standalone subproblem of the meta-problem
!# "FSI". In example 1 this is:
!#
!# t_problem("FSI")
!#   p_rproblemFirst => t_problem("Fluid")
!#   p_rproblemLast  => t_problem("Structure")
!#
!# These links make it possible to iterate over all subproblems given
!# the meta-problem as top-level problem structure.
!#
!# So far, the two subproblems are not releated to each other. This is
!# were p_rproblemPrev/Next comes into play. The relation between thw
!# the standalone subproblems "Fluid" and "Structure" is set up by the
!# initialisation procedure.
!#
!# t_problem("Fluid")
!#   p_rproblemTop  => t_problem("FSI")
!#   p_rproblemPrev => null()
!#   p_rproblemNext => t_problem("Structure")
!#
!# t_problem("Structure")
!#   p_rproblemTop  => t_problem("FSI")
!#   p_rproblemPrev => t_problem("Fluid")
!#   p_rproblemNext => null()
!#
!# These link make it possible to quickly change from within (!) one
!# standalone subproblem to the other (flat/horizontal hierarchy) and
!# to go to the top-level FSI problem (vertical hierarchy).
!#
!# Following example 2, the standalone structure subproblem has a
!# boundary parametrisation (in 2d) and a hierarchy of problem levels
!# which store matrices, vectors, discretisation structures, and the
!# triangulation on the particular refinement level. This again
!# vertical hierarchy is addressed, e.g., as follows:
!#
!# t_problem("Structure")
!#   p_rproblemLevelMin => t_problemLevel("Structure_level1") 
!#   p_rproblemLevelMax => t_problemLevel("Structure_level4")
!#
!# Each level knows its top and bottom neighbor and has a link to the
!# problem structure, e.g.
!#
!# t_problemLevel("Structure_level1")
!#   p_rproblemLevelCoarse => null()
!#   p_rproblemLevelFine   => t_problemLevel("Structure_level2")
!#   p_rproblem            => t_problem("Structure")
!#
!# t_problemLevel("Structure_level2")
!#   p_rproblemLevelCoarse => t_problemLevel("Structure_level1")
!#   p_rproblemLevelFine   => t_problemLevel("Structure_level3")
!#   p_rproblem            => t_problem("Structure")
!#
!# ...
!#
!# t_problemLevel("Structure_level4")
!#   p_rproblemLevelCoarse => t_problemLevel("Structure_level3")
!#   p_rproblemLevelFine   => null()
!#   p_rproblem            => t_problem("Structure")
!#
!# 
!# A similar substructure is attached to t_problem("Fluid"). However,
!# since the fluid problem consists of three embedded subproblems, the
!# situation is slightly different:
!#
!# t_problem("Fluid")
!#   p_rproblemFirst    => t_problem("compEuler1")
!#   p_rproblemLast     => t_problem("Transport")
!#   p_rproblemLevelMin => t_problemLevel("Fluid_level1") 
!#   p_rproblemLevelMax => t_problemLevel("Fluid_level4")
!#
!#   t_problem("compEuler1")
!#     p_rproblemTop      => t_problem("Fluid")
!#     p_rproblemPrev     => null()
!#     p_rproblemNext     => t_problem("compEuler2")
!#     p_rproblemLevelMin => t_problemLevel("compEuler1_level1") 
!#     p_rproblemLevelMax => t_problemLevel("compEuler1_level4")
!#
!#   t_problem("compEuler2")
!#     p_rproblemTop      => t_problem("Fluid")
!#     p_rproblemPrev     => t_problem("compEuler1")
!#     p_rproblemNext     => t_problem("Transport")
!#     p_rproblemLevelMin => t_problemLevel("compEuler2_level1") 
!#     p_rproblemLevelMax => t_problemLevel("compEuler2_level4")
!#
!#   t_problem("Transport")
!#     p_rproblemTop      => t_problem("Fluid")
!#     p_rproblemPrev     => t_problem("compEuler2")
!#     p_rproblemNext     => null()
!#     p_rproblemLevelMin => t_problemLevel("Transport_level1") 
!#     p_rproblemLevelMax => t_problemLevel("Transport_level4")
!#
!# Each embedded subproblem has its own level hierarchy and can work
!# on it individually. However, user making use of embedded
!# subproblems must be aware of the fact that some content is shared
!# (e.g., triangulation, boundary parametrisation) and that
!# modifications in one embedded subproblem affect all other
!# subproblems as well. This approach is quite fragile since removing
!# a single link may suffice to destroy the entire problem but it
!# admits utmost flexibility and code reuse. Declaring all
!# substructures private would prevent breaking the entire structure
!# but then each access to the content would require using some
!# get/set routine which is quite tedious.
!#
!#
!# The following routines are available:
!#
!# 1.) problem_initProblem
!#     -> Initialises a new problem structure
!#
!# 2.) problem_createProblem
!#     -> Creates a new problem structure
!#
!# 3.) problem_releaseProblem
!#     -> Releases an existing problem structure
!#
!# 4.) problem_appendProblem
!#     -> Appends a problem structure to the linked list of problem structures
!#
!# 5.) problem_prependProblem
!#     -> Prepends a problem structure to the linked list of problem structures
!#
!# 6.) problem_removeProblem
!#     -> Removes a problem structure from the linked list of problem structures
!#
!# 7.) problem_infoProblem / problem_infoProblemDescriptor
!#     -> Outputs information about the problem structure or its descriptor
!#
!# 8.) problem_createLevel
!#     -> Creates a new problem level structure
!#
!# 9.) problem_releaseLevel
!#     -> Releases an existing problem level structure
!#
!# 10.) problem_appendLevel
!#     -> Appends a problem level structure into the linked list of problem
!#        level structures
!#
!# 11.) problem_prependLevel
!#     -> Prepends a problem level structure into the linked list of problem
!#        level structures
!#
!# 12.) problem_removeLevel
!#      -> Removes an existing problem level structure from the linked list of
!#         problem level structures
!#
!# 13.) problem_infoLevel
!#      -> Outputs information about the problem level structure
!#
!# 14.) problem_getLevel = problem_getLevelDirect /
!#                         problem_getLevelIndirect
!#      -> Returns a pointer to the problem level specified by the level number
!#
!# 15.) problem_setSpec
!#      -> Sets the problem specification of one or more levels
!#
!# 16.) problem_initDescriptor
!#      -> Initialises a problem descriptor
!#
!# 17.) problem_combineDescriptors
!#      -> Combine two different problem descriptors
!#
!# 18.) problem_releaseDescriptor
!#      -> Releases a problem descriptor
!#
!#
!# The following auxiliary routines are available
!#
!# 1.) problem_initDiscrAll
!#     -> Initialises all discretisation structures of the problem
!#
!# 2.) problem_updateDiscrAll
!#     -> Updates all discretisation structures of the problem
!#
!# 3.) problem_initDiscr
!#     -> Initialises a single discretisation structure
!#
!# 4.) problem_updateDiscr
!#     -> Updates a single discretisation structures
!#
!# 5.) problem_initCubInfoAll
!#     -> Initialises all cubature info structures of the problem
!#
!# 6.) problem_updateCubInfoAll
!#     -> Updates all cubature info structures of the problem
!#
!# 7.) problem_initCubInfo
!#     -> Initialises a single cubature info structures
!#
!# 8.) problem_updateCubInfo
!#     -> Updates a single cubature info structures
!#
!# 9.) problem_initMatrixAll
!#     -> Initialises all scalar/block matrices of the problem
!#
!# 10.) problem_updateMatrixAll
!#     -> Updates all scalar/block matrices of the problem
!#
!# 11.) problem_initMatrixScalar / problem_initMatrixBlock
!#     -> Initialises a single scalar/block matrix
!#
!# 12.) problem_updateMatrixScalar / problem_updateMatrixBlock
!#      -> Updates a single scalar/block matrix
!#
!# 13.) problem_updateMatrix
!#      -> Updates a single scalar/block matrix
!#
!# 14.) problem_initVectorAll
!#     -> Initialises all scalar/block vectors of the problem
!#
!# 15.) problem_updateVectorAll
!#     -> Updates all scalar/block vectors of the problem
!#
!# 16.) problem_initVectorScalar / problem_initVectorBlock
!#     -> Initialises a single scalar/block vector
!#
!# 17.) problem_updateVectorScalar / problem_updateVecBl
!#      -> Updates a single scalar/block vector
!#
!# 18.) problem_updateVector
!#      -> Updates a single scalar/block vector
!#
!# 19.) problem_initAFCstabAll
!#     -> Initialises all stabilisation structures of AFC-type of the problem
!#
!# 20.) problem_updateAFCstabAll
!#     -> Updates all stabilisation structures of AFC-type of the problem
!#
!# 21.) problem_initAFCstab
!#     -> Initialises a single stabilisation structure of AFC-type
!#
!# 22.) problem_updateAFCstab
!#      -> Updates a single stabilisation structure of AFC-type
!#
!# 23.) problem_initGroupFEMBlockAll
!#     -> Initialises all blocks of group finite element structures of
!#        the problem
!#
!# 24.) problem_updateGroupFEMBlockAll
!#     -> Updates all blocks of group finite element structures of the
!#        problem
!#
!# 25.) problem_initGroupFEMBlock
!#     -> Initialises a single block of group finite element structures
!#
!# 26.) problem_updateGroupFEMBlock
!#      -> Updates a single block of group finite element structures
!#
!# 27.) problem_clearTaskList
!#     -> Releases a task list
!#
!# 28.) problem_infoTaskList
!#      -> Outputs information about the task list
!#
!# </purpose>
!##############################################################################

module problem

#include "../flagship.h"

!$use omp_lib
  use afcstabbase
  use basicgeometry
  use bilinearformevaluation
  use boundary
  use boundaryfilter
  use collection
  use cubature
  use derivatives
  use element
  use fparser
  use fsystem
  use genoutput
  use groupfembase
  use io
  use linearformevaluation
  use lineariser
  use linearsystemblock
  use linearsystemscalar
  use meshmodification
  use paramlist
  use scalarpde
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none

  private
  public :: t_problem
  public :: t_problemLevel
  public :: t_problemDescriptor
  public :: problem_initProblem
  public :: problem_createProblem
  public :: problem_releaseProblem
  public :: problem_appendProblem
  public :: problem_prependProblem
  public :: problem_removeProblem
  public :: problem_infoProblem
  public :: problem_infoProblemDescriptor
  public :: problem_createLevel
  public :: problem_releaseLevel
  public :: problem_appendLevel
  public :: problem_prependLevel
  public :: problem_removeLevel
  public :: problem_infoLevel
  public :: problem_getLevel
  public :: problem_setSpec
  public :: problem_initDescriptor
  public :: problem_combineDescriptors
  public :: problem_releaseDescriptor

  interface problem_infoTaskList
    module procedure problem_infoDiscrTL
    module procedure problem_infoCubInfoTL
    module procedure problem_infoMatrixTL
    module procedure problem_infoVectorTL
    module procedure problem_infoAFCstabTL
    module procedure problem_infoGroupFEMBlockTL
  end interface problem_infoTaskList

  public :: problem_infoTaskList

  interface problem_clearTaskList
    module procedure problem_clearDiscrTL
    module procedure problem_clearCubInfoTL
    module procedure problem_clearMatrixTL
    module procedure problem_clearVectorTL
    module procedure problem_clearAFCstabTL
    module procedure problem_clearGroupFEMBlockTL
  end interface problem_clearTaskList

  public :: problem_clearTaskList

  public :: t_discrTask
  public :: problem_initDiscrAll
  public :: problem_updateDiscrAll
  public :: problem_initDiscr
  public :: problem_updateDiscr

  interface problem_updateCubInfo
    module procedure problem_updateCubInfo1
    module procedure problem_updateCubInfo2
  end interface problem_updateCubInfo

  public :: t_cubinfoTask
  public :: problem_initCubInfoAll
  public :: problem_updateCubInfoAll
  public :: problem_initCubInfo
  public :: problem_updateCubInfo
  
  public :: t_matrixTask
  public :: problem_initMatrixAll
  public :: problem_updateMatrixAll
  public :: problem_initMatrixScalar
  public :: problem_initMatrixBlock
  public :: problem_updateMatrixScalar
  public :: problem_updateMatrixBlock
  public :: problem_updateMatrix

  public :: t_vectorTask
  public :: problem_initVectorAll
  public :: problem_updateVectorAll
  public :: problem_initVectorScalar
  public :: problem_initVectorBlock
  public :: problem_updateVectorScalar
  public :: problem_updateVecBl
  public :: problem_updateVector

  interface problem_updateAFCstab
    module procedure problem_updateAFCstab1
    module procedure problem_updateAFCstab2
  end interface problem_updateAFCstab

  public :: t_afcstabTask
  public :: problem_initAFCstabAll
  public :: problem_updateAFCstabAll
  public :: problem_initAFCstab
  public :: problem_updateAFCstab

  interface problem_updateGroupFEMBlock
    module procedure problem_updateGroupFEMBlock1
    module procedure problem_updateGroupFEMBlock2
  end interface problem_updateGroupFEMBlock

  public :: t_groupFEMBlockTask
  public :: problem_initGroupFEMBlockAll
  public :: problem_updateGroupFEMBlockAll
  public :: problem_initGroupFEMBlock
  public :: problem_updateGroupFEMBlock

  !*****************************************************************************

  interface problem_initProblem
    module procedure problem_initProblemDirect
    module procedure problem_initProblemFromParlst
  end interface problem_initProblem

  interface problem_getLevel
    module procedure problem_getLevelDirect
    module procedure problem_getLevelIndirect
  end interface problem_getLevel

  !*****************************************************************************

!<constants>

!<constantblock description="Flags to specify the problem type">

  ! Unspecified problem
  integer(I32), parameter, public :: PROBTYPE_UNDEFINED              = 0

  ! Meta-problem (this means that it does not provide
  ! data itself is will most likely not work as top-level
  ! problem within a standard application)
  integer(I32), parameter, public :: PROBTYPE_META                   = 2**0

  ! Standalone problem (this means that it can be
  ! used as top-level problem within an application)
  integer(I32), parameter, public :: PROBTYPE_STANDALONE             = 2**1

  ! Embedded problem (this means that it is part of a larger
  ! problem; its content is shared with the top-level problem)
  integer(I32), parameter, public :: PROBTYPE_EMBEDDED               = 2**2

!</constantblock>

!<constantblock description="Flags for the problem descriptor specification bitfield">

  ! Standard problem
  integer(I32), parameter, public :: PROBDESC_PSPEC_STANDARD         = 0

  ! Generate extended raw mesh
  integer(I32), parameter, public :: PROBDESC_PSPEC_EXTENDEDRAW      = 2**0

!</constantblock>

!<constantblock description="Flags for the action">

  ! Action: none
  integer, parameter, public :: PROBACTION_NONE                      = 0

  ! Action: create
  integer, parameter, public :: PROBACTION_CREATE                    = 1

  ! Action: duplicate
  integer, parameter, public :: PROBACTION_DUPLICATE                 = 2

!</constantblock>

!<constantblock description="Flags for the action performance">

  ! Action perform: never
  integer(I32), parameter, public :: PROBACTION_PERFORM_NEVER        = 0

  ! Action perform: on init
  integer(I32), parameter, public :: PROBACTION_PERFORM_INIT         = 2**0

  ! Action perform: on update
  integer(I32), parameter, public :: PROBACTION_PERFORM_UPDATE       = 2**1

  ! Action perform: always
  integer(I32), parameter, public :: PROBACTION_PERFORM_ALWAYS       = PROBACTION_PERFORM_INIT&
                                                                     + PROBACTION_PERFORM_UPDATE

!</constantblock>

!<constantblock description="Flags for the method to be used for matrix creation">

  ! Create virtual matrix without data
  integer, parameter,  public :: PROBMETH_MATRIX_VIRTUAL             = 0

  ! Create empty matrix
  integer, parameter,  public :: PROBMETH_MATRIX_EMPTY               = 1

  ! Create identify matrix
  integer, parameter,  public :: PROBMETH_MATRIX_IDENTITY            = 2
  
  ! Create matrix from bilinear form
  integer, parameter,  public :: PROBMETH_MATRIX_BILF                = 3

  ! Create matrix by standard lumping of another matrix
  integer, parameter,  public :: PROBMETH_MATRIX_LUMP_STD            = 4

  ! Create matrix by diagonal lumping of another matrix
  integer, parameter,  public :: PROBMETH_MATRIX_LUMP_DIAG           = 5

  ! Create matrix with extended sparsity pattern
  integer, parameter,  public :: PROBMETH_MATRIX_EXTENDSPARSITY      = 6

!</constantblock>

!<constantblock description="Flags for the method to be used for vector creation">

  ! Create virtual vector without data
  integer, parameter,  public :: PROBMETH_VECTOR_VIRTUAL             = 0

  ! Create empty vector
  integer, parameter,  public :: PROBMETH_VECTOR_EMPTY               = 1

  ! Create unit vector
  integer, parameter,  public :: PROBMETH_VECTOR_UNITY               = 2
  
  ! Create vector from bilinear form
  integer, parameter,  public :: PROBMETH_VECTOR_LINF                = 3
  
  ! Create coordinates of the degrees of freedom
  integer, parameter,  public :: PROBMETH_VECTOR_DOF                 = 4

!</constantblock>

!</constants>

  !*****************************************************************************

!<types>

!<typeblock>

  ! This data structure contains the abstract problem description

  type t_problemDescriptor

    ! Problem descriptor specification tag. This is a bitfield coming
    ! from an OR combination of different PROBDESC_PSPEC_xxxx
    ! constants and specifies various details of the problem
    ! descriptor. If it is =PROBDESC_PSPEC_STANDARD, the problem
    ! descriptor is a usual descriptor that needs no special handling.
    integer(I32) :: iproblemSpec = PROBDESC_PSPEC_STANDARD

    ! --------------------------------------------------------------------------
    ! The following data must be provided for all problem types

    ! Name of the problem
    character(len=SYS_STRLEN) :: cproblem = ''
    
    ! Type of the problem
    integer(I32) :: cproblemtype = PROBTYPE_UNDEFINED

    ! If this descriptor represents a meta-problem which does not
    ! provide own content then we must have an array of subproblems.
    ! If this descriptor represents a standalone problem with one or
    ! more embedded subproblems then we must have an array of
    ! subproblems, too. The difference between meta-problems and
    ! standalone problems is that meta-problems may not have any data,
    ! whereas standalone problems with possible subproblems do have
    ! their own problem level structure.
    type(t_problemDescriptor), dimension(:), pointer :: Rsubproblem => null()
    
    ! --------------------------------------------------------------------------
    ! STANDALONE PROBLEM DESCRIPTOR
    !
    ! If this descriptor represents a standalone problem the following
    ! data must be provided by the user to generate the level structure.

    ! Spatial dimension
    integer :: ndimension = 0

    ! Minimum problem level
    integer :: nlmin = 0

    ! Maximum problem level
    integer :: nlmax = 0

    ! Number of discretisation structures
    integer :: ndiscretisation = 0

    ! Number of cubature info structures
    integer :: ncubatureinfo = 0

    ! Number of AFC stabilisations
    integer :: nafcstab = 0

    ! Number of group finite element blocks
    integer :: ngroupfemBlock = 0

    ! Number of scalar matrices
    integer :: nmatrixScalar = 0

    ! Number of block matrices
    integer :: nmatrixBlock = 0

    ! Number of scalar vectors
    integer :: nvectorScalar = 0

    ! Number of block vectors
    integer :: nvectorBlock = 0

#ifdef ENABLE_COPROCESSOR_SUPPORT
    ! Number of CUDA streams
    integer :: nstream = 0
#endif

    ! Name of the triangulation file
    character(LEN=SYS_STRLEN) :: trifile = ''

    ! Name of the boundary parametrisation file
    character(LEN=SYS_STRLEN) :: prmfile = ''

    ! Strategy to convert triangulation
    integer :: iconvStrategy = TRI_CONVERT_NONE

    ! Amount of stochastically mesh disturbance
    real(DP) :: ddisturbMesh = 0.0_DP

  end type t_problemDescriptor

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure contains the top-level structure of a complete
  ! problem configuration. In a multi-physics problem, this structure
  ! serves as meta-problem which contains multiple subproblems. Each
  ! subproblem can itself be a meta-problem containing multiple subproblems.
  !
  ! Each real problem may contain a single problem level or a sequence of
  ! nested problem levels required for applying multigrid-type solvers.

  type t_problem

    ! Name of the problem
    character(len=SYS_STRLEN) :: cproblem = ''

    ! Type of the problem
    integer(I32) :: cproblemtype = PROBTYPE_UNDEFINED

    ! If this structure is used to represent a meta-problem
    ! then only the two following quantities apply.
    ! --------------------------------------------------------------------------

    ! Pointer to the first subproblem instance which is part of this problem
    type(t_problem), pointer :: p_rproblemFirst => null()

    ! Pointer to the last subproblem instance which is part of this problem
    type(t_problem), pointer :: p_rproblemLast => null()

    ! If this structure is used to represent a real problem
    ! then the following quantities apply.
    ! --------------------------------------------------------------------------

    ! Pointer to the outer top-level problem
    type(t_problem), pointer :: p_rproblemTop => null()

    ! Pointer to the previous problem instance
    type(t_problem), pointer :: p_rproblemPrev => null()

    ! Pointer to the next problem instance
    type(t_problem), pointer :: p_rproblemNext => null()

    ! Boundary parametrisation which is the same for all problem
    ! levels. A boundary parametrisation exists in 2D only and will
    ! not contain any data in 1D and 3D.
    type(t_boundary), pointer :: rboundary => null()

    ! Pointers to the minimum problem level
    type(t_problemLevel), pointer :: p_rproblemLevelMin => null()

    ! Pointers to the maximum problem level
    type(t_problemLevel), pointer :: p_rproblemLevelMax => null()

    ! --------------------------------------------------------------------------

    ! INTERNAL DATA: so-called task lists can be attached to the
    ! problem structure which represent the dependency between tasks
    ! such as create, share, copy of the content
    type(t_discrTask), pointer :: p_rdiscrTasklist => null()

    ! INTERNAL DATA: so-called task lists can be attached to the
    ! problem structure which represent the dependency between tasks
    ! such as create, share, copy of the content
    type(t_cubinfoTask), pointer :: p_rcubinfoTasklist => null()

    ! INTERNAL DATA: so-called task lists can be attached to the
    ! problem structure which represent the dependency between tasks
    ! such as create, share, copy of the content
    type(t_matrixTask), pointer :: p_rmatrixTasklist => null()

    ! INTERNAL DATA: so-called task lists can be attached to the
    ! problem structure which represent the dependency between tasks
    ! such as create, share, copy of the content
    type(t_vectorTask), pointer :: p_rvectorTasklist => null()

    ! INTERNAL DATA: so-called task lists can be attached to the
    ! problem structure which represent the dependency between tasks
    ! such as create, share, copy of the content
    type(t_afcstabTask), pointer :: p_rafcstabTasklist => null()

    ! INTERNAL DATA: so-called task lists can be attached to the
    ! problem structure which represent the dependency between tasks
    ! such as create, share, copy of the content
    type(t_groupFEMBlockTask), pointer :: p_rgroupFEMBlockTasklist => null()

  end type t_problem

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure contains all required data on one problem level

  type t_problemLevel

    ! Specification tag. This bitfield can be used by the application
    ! programmer to specify various details of the problem level.
    integer(I32) :: iproblemSpec = 0_I32

    ! Number of the problem level
    integer :: ilev

#ifdef ENABLE_COPROCESSOR_SUPPORT
    ! Array of CUDA streams
    integer(I64), dimension(:), pointer :: Istream => null()
#endif

    ! Triangulation structure
    type(t_triangulation), pointer :: rtriangulation

    ! Array of discretisation structure
    type(t_blockDiscretisation), dimension(:), pointer :: Rdiscretisation => null()

    ! Array of scalar cubature info structures
    type(t_scalarCubatureInfo), dimension(:), pointer :: RcubatureInfo => null()

    ! Array of AFC stabilisations
    type(t_afcstab), dimension(:), pointer :: Rafcstab => null()

    ! Array of group finite element blocks
    type(t_groupFEMBlock), dimension(:), pointer :: RgroupFEMBlock => null()

    ! Array of scalar matrices
    type(t_matrixScalar), dimension(:), pointer :: RmatrixScalar => null()

    ! Array of block matrices
    type(t_matrixBlock), dimension(:), pointer :: RmatrixBlock => null()

    ! Array of scalar vectors
    type(t_vectorScalar), dimension(:), pointer :: RvectorScalar => null()

    ! Array of block vectors
    type(t_vectorBlock), dimension(:), pointer :: RvectorBlock => null()

    ! Pointer to the encompassing problem structure
    type(t_problem), pointer :: p_rproblem => null()

    ! Pointer to next coarse problem level
    type(t_problemLevel), pointer :: p_rproblemLevelCoarse => null()

    ! Pointer to next finer problem level
    type(t_problemLevel), pointer :: p_rproblemLevelFine => null()

  end type t_problemLevel

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure represents a single task item for block discretisations
  type t_discrTask

    ! Task to be performed: CREATE, DUPLICATE
    integer :: ctask = PROBACTION_NONE

    ! Task specification tag. This is a bitfield coming from an OR
    ! combination of different PROBACTION_PERFORM_xxxx constants and
    ! specifies when to perform the particular action.
    integer(i32) :: iperform = PROBACTION_PERFORM_NEVER

    ! Flag to define wether structure should be copied or shared
    logical :: bshareStructure = .true.

    ! First block in the discretisation structure to duplicate
    integer :: ifirstBlock = 0

    ! Last block in the discretisation structure to duplicate
    integer :: ilastBlock = 0
    
    ! Name of the section in the parameter list which provides all
    ! necessary information to create/update this task
    character(len=SYS_STRLEN) :: ssectionName = ''

    ! Pointer to the destination problem
    type(t_problem), pointer :: p_rproblemDest => null()

    ! Pointer to the source problem
    type(t_problem), pointer :: p_rproblemSrc => null()

    ! Pointer to the destination problem level
    type(t_problemLevel), pointer :: p_rproblemLevelDest => null()

    ! Pointer to the sourceproblem level
    type(t_problemLevel), pointer :: p_rproblemLevelSrc => null()

    ! Pointer to the destination block discretisation structure
    type(t_blockDiscretisation), pointer :: p_rblockDiscrDest => null()

    ! Pointer to the source block discretisation structure (if any)
    type(t_blockDiscretisation), pointer :: p_rblockDiscrSrc  => null()

    ! Pointer to the next task to be performed after this task
    type(t_discrTask), pointer :: p_rnextTask => null()

  end type t_discrTask
  
!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure represents a single task item for cubature info structures
  type t_cubinfoTask

    ! Task to be performed: CREATE, DUPLICATE
    integer :: ctask = PROBACTION_NONE
    
    ! Task specification tag. This is a bitfield coming from an OR
    ! combination of different PROBACTION_PERFORM_xxxx constants and
    ! specifies when to perform the particular action.
    integer(i32) :: iperform = PROBACTION_PERFORM_NEVER
    
    ! Name of the section in the parameter list which provides all
    ! necessary information to create/update this task
    character(len=SYS_STRLEN) :: ssectionName = ''

    ! Pointer to the destination problem
    type(t_problem), pointer :: p_rproblemDest => null()

    ! Pointer to the source problem
    type(t_problem), pointer :: p_rproblemSrc => null()

    ! Pointer to the destination problem level
    type(t_problemLevel), pointer :: p_rproblemLevelDest => null()

    ! Pointer to the sourceproblem level
    type(t_problemLevel), pointer :: p_rproblemLevelSrc => null()

    ! Pointer to the destination cubature info structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfoDest => null()

    ! Pointer to the source cubature info structure (if any)
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfoSrc  => null()

    ! Pointer to the next task to be performed after this task
    type(t_cubinfoTask), pointer :: p_rnextTask => null()

    !---------- INTERNAL DATA --------------------------------------------------
    
    ! Type of cubature rule
    integer(I32) :: ccubType = CUB_UNDEFINED
    
    ! Number of levels for summed cubature
    integer :: nlevels = 0
    
    ! Pointer to the spatial discretisation structure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr => null()

  end type t_cubinfoTask

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure represents a single task item for matrices
  type t_matrixTask

    ! Flag to identify scalar and block matrices
    logical :: bisMatrixScalar = .true.

    ! Task to be performed: CREATE, DUPLICATE
    integer :: ctask = PROBACTION_NONE
    
    ! Task specification tag. This is a bitfield coming from an OR
    ! combination of different PROBACTION_PERFORM_xxxx constants and
    ! specifies when to perform the particular action.
    integer(i32) :: iperform = PROBACTION_PERFORM_NEVER

    ! Flag to define wether structure should be copied or shared
    integer :: cdupStructure = LSYSSC_DUP_IGNORE

    ! Flag to define wether content should be copied or shared
    integer :: cdupContent = LSYSSC_DUP_IGNORE

    ! Name of the section in the parameter list which provides all
    ! necessary information to create/update this task
    character(len=SYS_STRLEN) :: ssectionName = ''

    ! Pointer to the destination problem
    type(t_problem), pointer :: p_rproblemDest => null()

    ! Pointer to the source problem
    type(t_problem), pointer :: p_rproblemSrc => null()

    ! Pointer to the destination problem level
    type(t_problemLevel), pointer :: p_rproblemLevelDest => null()

    ! Pointer to the sourceproblem level
    type(t_problemLevel), pointer :: p_rproblemLevelSrc => null()

    ! Pointer to the destination scalar matrix
    type(t_matrixScalar), pointer :: p_rmatrixScalarDest => null()

    ! Pointer to the source scalar matrix (if any)
    type(t_matrixScalar), pointer :: p_rmatrixScalarSrc  => null()
    
    ! Pointer to the destination block matrix
    type(t_matrixBlock), pointer :: p_rmatrixBlockDest => null()

    ! Pointer to the source block matrix (if any)
    type(t_matrixBlock), pointer :: p_rmatrixBlockSrc  => null()

    ! Pointer to the next task to be performed after this task
    type(t_matrixTask), pointer :: p_rnextTask => null()

    !---------- INTERNAL DATA --------------------------------------------------

    ! Method to be used: PROBMETH_MATRIX_xxx
    integer, dimension(:), pointer :: Cmethod => null()

    !---------- INTERNAL DATA FOR SCALAR MATRIX --------------------------------

    ! Data type of the entries in the matrix
    integer :: cdataType = ST_NOHANDLE

    ! Matrix format can take one of the LSYSSC_MATRIXx format flags
    integer :: cmatrixFormat = LSYSSC_MATRIXUNDEFINED

    ! Interleaved matrix format can take LSYSSC_MATRIX1 or LSYSSC_MATRIXD
    integer :: cinterleavematrixFormat = LSYSSC_MATRIXUNDEFINED

    ! Number of variables for interleaved matrix
    integer :: nvar = 0

    ! Scaling factor
    real(DP) :: dscaleFactor = 0.0_DP

    ! Pointer to the spatial discretisation for the trial functions.
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrTrial => null()

    ! Pointer to the spatial discretisation for the test functions.
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrTest => null()

    ! Pointer to the cubature info structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo => null()

    ! Pointer to the scalar bilinear form 
    type(t_bilinearForm), pointer :: p_rform => null()

    ! Pointer to a function parser
    type(t_fparser), pointer :: p_rfparser => null()
    
    !---------- INTERNAL DATA FOR BLOCK MATRIX ---------------------------------

    ! Pointer to the block discretisation for the trial functions.
    type(t_blockDiscretisation), pointer :: p_rblockDiscrTrial => null()

    ! Pointer to the block discretisation for the test functions.
    type(t_blockDiscretisation), pointer :: p_rblockDiscrTest => null()

  end type t_matrixTask

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure represents a single task item for vectors
  type t_vectorTask

    ! Flag to identify scalar and block vectors
    logical :: bisVectorScalar = .true.

    ! Task to be performed: CREATE, DUPLICATE
    integer :: ctask = PROBACTION_NONE
    
    ! Task specification tag. This is a bitfield coming from an OR
    ! combination of different PROBACTION_PERFORM_xxxx constants and
    ! specifies when to perform the particular action.
    integer(i32) :: iperform = PROBACTION_PERFORM_NEVER

    ! Flag to define wether structure should be copied or shared
    integer :: cdupStructure = LSYSSC_DUP_IGNORE

    ! Flag to define wether content should be copied or shared
    integer :: cdupContent = LSYSSC_DUP_IGNORE

    ! Name of the section in the parameter list which provides all
    ! necessary information to create/update this task
    character(len=SYS_STRLEN) :: ssectionName = ''

    ! Pointer to the destination problem
    type(t_problem), pointer :: p_rproblemDest => null()

    ! Pointer to the source problem
    type(t_problem), pointer :: p_rproblemSrc => null()

    ! Pointer to the destination problem level
    type(t_problemLevel), pointer :: p_rproblemLevelDest => null()

    ! Pointer to the sourceproblem level
    type(t_problemLevel), pointer :: p_rproblemLevelSrc => null()

    ! Pointer to the destination scalar vector
    type(t_vectorScalar), pointer :: p_rvectorScalarDest => null()

    ! Pointer to the source scalar vector (if any)
    type(t_vectorScalar), pointer :: p_rvectorScalarSrc  => null()
    
    ! Pointer to the destination block vector
    type(t_vectorBlock), pointer :: p_rvectorBlockDest => null()

    ! Pointer to the source block vector (if any)
    type(t_vectorBlock), pointer :: p_rvectorBlockSrc  => null()

    ! Pointer to the next task to be performed after this task
    type(t_vectorTask), pointer :: p_rnextTask => null()

    !---------- INTERNAL DATA --------------------------------------------------

    ! Method to be used: PROBMETH_VECTOR_xxx
    integer, dimension(:), pointer :: Cmethod => null()

    !---------- INTERNAL DATA FOR SCALAR VECTOR --------------------------------

    ! Data type of the entries in the vector
    integer :: cdataType = ST_NOHANDLE

    ! Number of variables for interleaved vector
    integer :: nvar = 0

    ! Pointer to the spatial discretisation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr => null()

    ! Pointer to the cubature info structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo => null()

    ! Pointer to the scalar linear form 
    type(t_linearForm), pointer :: p_rform => null()

    ! Pointer to a function parser
    type(t_fparser), pointer :: p_rfparser => null()
    
    !---------- INTERNAL DATA FOR BLOCK VECTOR ---------------------------------

    ! Pointer to the block discretisation
    type(t_blockDiscretisation), pointer :: p_rblockDiscr => null()

  end type t_vectorTask

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure represents a single task item for
  ! stabilisation structures of AFC-type
  type t_afcstabTask

    ! Task to be performed: CREATE, DUPLICATE
    integer :: ctask = PROBACTION_NONE
    
    ! Task specification tag. This is a bitfield coming from an OR
    ! combination of different PROBACTION_PERFORM_xxxx constants and
    ! specifies when to perform the particular action.
    integer(i32) :: iperform = PROBACTION_PERFORM_NEVER
    
    ! Name of the section in the parameter list which provides all
    ! necessary information to create/update this task
    character(len=SYS_STRLEN) :: ssectionName = ''

    ! Pointer to the destination problem
    type(t_problem), pointer :: p_rproblemDest => null()

    ! Pointer to the source problem
    type(t_problem), pointer :: p_rproblemSrc => null()

    ! Pointer to the destination problem level
    type(t_problemLevel), pointer :: p_rproblemLevelDest => null()

    ! Pointer to the source problem level
    type(t_problemLevel), pointer :: p_rproblemLevelSrc => null()

    ! Pointer to the destination stabilisation structure
    type(t_afcstab), pointer :: p_rafcstabDest => null()

    ! Pointer to the source stabilisation structure (if any)
    type(t_afcstab), pointer :: p_rafcstabSrc  => null()

    ! Pointer to the next task to be performed after this task
    type(t_afcstabTask), pointer :: p_rnextTask => null()

    !---------- INTERNAL DATA --------------------------------------------------

  end type t_afcstabTask

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure represents a single task item for
  ! blocks of group finite element structures
  type t_groupFEMBlockTask

    ! Task to be performed: CREATE, DUPLICATE
    integer :: ctask = PROBACTION_NONE
    
    ! Task specification tag. This is a bitfield coming from an OR
    ! combination of different PROBACTION_PERFORM_xxxx constants and
    ! specifies when to perform the particular action.
    integer(i32) :: iperform = PROBACTION_PERFORM_NEVER
    
    ! Name of the section in the parameter list which provides all
    ! necessary information to create/update this task
    character(len=SYS_STRLEN) :: ssectionName = ''

    ! Pointer to the destination problem
    type(t_problem), pointer :: p_rproblemDest => null()

    ! Pointer to the source problem
    type(t_problem), pointer :: p_rproblemSrc => null()

    ! Pointer to the destination problem level
    type(t_problemLevel), pointer :: p_rproblemLevelDest => null()

    ! Pointer to the source problem level
    type(t_problemLevel), pointer :: p_rproblemLevelSrc => null()

    ! Pointer to the destination block of group finite element structures
    type(t_groupFEMBlock), pointer :: p_rgroupFEMBlockDest => null()

    ! Pointer to the source block of group finite element structures (if any)
    type(t_groupFEMBlock), pointer :: p_rgroupFEMBlockSrc  => null()

    ! Pointer to the next task to be performed after this task
    type(t_groupFEMBlockTask), pointer :: p_rnextTask => null()

    !---------- INTERNAL DATA --------------------------------------------------

  end type t_groupFEMBlockTask

!</typeblock>

!</types>

contains

  !*****************************************************************************

!<subroutine>

  subroutine problem_initProblemDirect(rproblem, rproblemDescriptor)

!<description>
    ! This subroutine initialises an abstract problem structure from
    ! the given problem descriptor structure
!</description>

!<input>
    ! abstract problem descriptor
    type(t_problemDescriptor), intent(in) :: rproblemDescriptor
!</input>

!<output>
    ! problem data structure
    type(t_problem), intent(out) :: rproblem
!</output>
!</subroutine>

    ! Initialise problem structure recursively.
    call problem_createProblem(rproblem)
    call recursiveInit(rproblem, rproblemDescriptor)

  contains

    ! Here, the real initialisation routine follows.

    !**************************************************************
    ! Internal routine which initialises the problem structure
    recursive subroutine recursiveInit(rproblem, rproblemDescriptor)

      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
      type(t_problem), intent(inout) :: rproblem
      
      ! local variables
      type(t_problem), pointer :: rsubproblem
      type(t_problemLevel), pointer :: rproblemLevel,p_rproblemLevel,p_rtopLevel
      integer :: i
      logical :: bextendedRaw
      
      ! Initialise global problem structure
      rproblem%cproblem     = rproblemDescriptor%cproblem
      rproblem%cproblemtype = rproblemDescriptor%cproblemtype

      ! Do we have to initialise content for this subproblem?
      if (iand(rproblemDescriptor%cproblemtype,PROBTYPE_META) .eq. 0) then

        !-----------------------------------------------------------------------
        ! Part 1: Initialise the coarsest problem level structure.
        nullify(rproblemLevel); allocate(rproblemLevel)
        call problem_createLevel(rproblemLevel, rproblemDescriptor%nlmin)
        bextendedRaw = (iand(rproblemDescriptor%iproblemSpec,&
                             PROBDESC_PSPEC_EXTENDEDRAW) .eq. 0)
        
        if (iand(rproblemDescriptor%cproblemtype, PROBTYPE_EMBEDDED) .eq. 0) then
          
          ! Allocate triangulation structure
          allocate(rproblemLevel%rtriangulation)

          ! Get triangulation from descriptor
          select case(rproblemDescriptor%ndimension)
          case (NDIM1D)
            ! Read coarse mesh from TRI-file
            call tria_readTriFile1D(rproblemLevel%rtriangulation,&
                rproblemDescriptor%trifile, bextendedRaw)

            ! Refine coarse mesh to minimum problem level
            call tria_quickRefine2LevelOrdering(rproblemDescriptor%nlmin-1,&
                rproblemLevel%rtriangulation)

          case (NDIM2D)
            ! Create new boundary and read from PRM-file
            allocate(rproblem%rboundary)
            call boundary_read_prm(rproblem%rboundary,&
                rproblemDescriptor%prmfile)

            ! Read coarse mesh from TRI-file
            call tria_readTriFile2D(rproblemLevel%rtriangulation,&
                rproblemDescriptor%trifile, rproblem%rboundary, .true.)

            ! Convert quadrilaterals to triangules if required
            if (rproblemDescriptor%iconvStrategy .ne. TRI_CONVERT_NONE) then
              call tria_rawGridToTri(rproblemLevel%rtriangulation,&
                  rproblemDescriptor%iconvStrategy)
            end if

            ! Generate extended raw mesh
            if (bextendedRaw) call tria_initExtendedRawMesh(rproblemLevel%rtriangulation)

            ! Refine coarse mesh to minimum problem level using an analytic
            ! boundary parametrisation
            call tria_quickRefine2LevelOrdering(rproblemDescriptor%nlmin-1,&
                rproblemLevel%rtriangulation, rproblem%rboundary)

          case (NDIM3D)
            ! Read coarse mesh from TRI-file
            call tria_readTriFile3D(rproblemLevel%rtriangulation,&
                rproblemDescriptor%trifile, bnoExtendedRaw=.true.)

            ! Convert hexahedrals to tetrahedrals if required
            if (rproblemDescriptor%iconvStrategy .ne. TRI_CONVERT_NONE) then
              call output_line('Conversion to tetrahedrals is not available yet!',&
                  OU_CLASS_WARNING,OU_MODE_STD,'problem_initProblem')
              call sys_halt()
            end if

            ! Generate extended raw mesh
            if (bextendedRaw) call tria_initExtendedRawMesh(rproblemLevel%rtriangulation)

            ! Refine coarse mesh to minimum problem level
            call tria_quickRefine2LevelOrdering(rproblemDescriptor%nlmin-1,&
                rproblemLevel%rtriangulation)

          case default
            call output_line('Invalid spatial dimension!',&
                OU_CLASS_WARNING,OU_MODE_STD,'problem_initProblem')
            call sys_halt()
          end select

          ! Disturb mesh stochastically if required
          if (rproblemDescriptor%ddisturbmesh > 0.0_DP) then
            call meshmod_disturbMesh(rproblemLevel%rtriangulation,&
                rproblemDescriptor%ddisturbMesh)
          end if

          ! Create standard mesh from raw mesh
          if (associated(rproblem%rboundary)) then
            call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
                rproblem%rboundary)
          else
            call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation)
          end if
        else

          ! Get corresponding level from top-level problem
          p_rtopLevel => problem_getLevel(rproblem%p_rproblemTop,&
              rproblemLevel%ilev, .true.)
          if (.not.associated(p_rtopLevel)) then
            call output_line('Internal error in problem structure!',&
                OU_CLASS_ERROR, OU_MODE_STD, 'problem_initProblem')
            call sys_halt()
          end if

          ! Duplicate triangulation structure from top-level problem
          rproblemLevel%rtriangulation => p_rtopLevel%rtriangulation
        end if

        ! Allocate auxiliary arrays
        call allocateStructure(rproblemLevel, rproblemDescriptor)

        ! Append coarsest problem level to global problem structure
        call problem_appendLevel(rproblem, rproblemLevel)
        p_rproblemLevel => rproblemLevel

        !-----------------------------------------------------------------------
        ! Part 2: Initialise the finer problem level structures.
        do i = rproblemDescriptor%nlmin+1, rproblemDescriptor%nlmax

          ! Initialise current problem level structure
          nullify(rproblemLevel); allocate(rproblemLevel)
          call problem_createLevel(rproblemLevel, i)

          if (iand(rproblemDescriptor%cproblemtype, PROBTYPE_EMBEDDED) .eq. 0) then

            ! Allocate triangulation structure
            allocate(rproblemLevel%rtriangulation)

            if (associated(rproblem%rboundary)) then
              ! Generate regularly refined mesh from coarser level
              call tria_refine2LevelOrdering(p_rproblemLevel%rtriangulation,&
                  rproblemLevel%rtriangulation, rproblem%rboundary)
              
              ! Create standard mesh from raw mesh
              call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
                  rproblem%rboundary)
            else
              ! Generate regularly refined mesh from coarser level
              call tria_refine2LevelOrdering(p_rproblemLevel%rtriangulation,&
                  rproblemLevel%rtriangulation)
              
              ! Create standard mesh from raw mesh
              call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation)
            end if
          else

            ! Get corresponding level from top-level problem
            p_rtopLevel => problem_getLevel(rproblem%p_rproblemTop,&
                rproblemLevel%ilev, .true.)
            if (.not.associated(p_rtopLevel)) then
              call output_line('Internal error in problem structure!',&
                  OU_CLASS_ERROR, OU_MODE_STD, 'problem_initProblem')
              call sys_halt()
            end if

            ! Duplicate triangulation structure from top-level problem
            rproblemLevel%rtriangulation => p_rtopLevel%rtriangulation
          end if

          ! Allocate auxiliary arrays
          call allocateStructure(rproblemLevel, rproblemDescriptor)     

          ! Append current problem level to global problem structre
          call problem_appendLevel(rproblem, rproblemLevel)

          ! Proceed to next level
          p_rproblemLevel => rproblemLevel
        end do

        !-----------------------------------------------------------------------
        ! Part 3: Compress triangulation structure
        if (iand(rproblemDescriptor%cproblemtype, PROBTYPE_EMBEDDED) .eq. 0) then       
          p_rproblemLevel => rproblem%p_rproblemLevelMax
          do while (associated(p_rproblemLevel) .and.&
              associated(p_rproblemLevel%p_rproblemLevelCoarse))
            call tria_compress2LevelOrdHierarchy(&
                p_rproblemLevel%rtriangulation,&
                p_rproblemLevel%p_rproblemLevelCoarse%rtriangulation)
            p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
          end do
        end if

      end if ! end of meta-problem-check

      !-------------------------------------------------------------------------
      ! Part 4: Do we have to initialise any subproblems?
      if (associated(rproblemDescriptor%Rsubproblem)) then
        
        ! Initialise each subproblem recursively
        do i = lbound(rproblemDescriptor%Rsubproblem,1),&
               ubound(rproblemDescriptor%Rsubproblem,1)
          
          ! Ok, let us create a new problem structure ...
          nullify(rsubproblem); allocate(rsubproblem)
          call problem_createProblem(rsubproblem)
          
          ! ... append it to the list of subproblems ...
          call problem_appendProblem(rproblem, rsubproblem)          

          ! and initialise it recursively.
          call recursiveInit(rsubproblem, rproblemDescriptor%Rsubproblem(i))
        end do
        
        ! If this is a meta-problem we are done
        if (iand(rproblemDescriptor%cproblemtype,PROBTYPE_META) .ne. 0) return
      end if

    end subroutine recursiveInit

    !**************************************************************
    ! Internal routine which allocates all problem level structures
    subroutine allocateStructure(rproblemLevel, rproblemDescriptor)

      type(t_problemLevel), intent(inout) :: rproblemLevel
      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
#ifdef ENABLE_COPROCESSOR_SUPPORT
      integer :: istream
#endif

      ! Allocate matrices, vectors, stabilisations, etc.
      if (rproblemDescriptor%ndiscretisation .gt. 0)&
          allocate(rproblemLevel%Rdiscretisation(&
          rproblemDescriptor%ndiscretisation))

      if (rproblemDescriptor%ncubatureInfo .gt. 0)&
          allocate(rproblemLevel%RcubatureInfo(&
          rproblemDescriptor%ncubatureInfo))

      if (rproblemDescriptor%nafcstab .gt. 0)&
          allocate(rproblemLevel%Rafcstab(&
          rproblemDescriptor%nafcstab))

      if (rproblemDescriptor%ngroupfemBlock .gt. 0)&
          allocate(rproblemLevel%RgroupFEMBlock(&
          rproblemDescriptor%ngroupfemBlock))

      if (rproblemDescriptor%nmatrixScalar .gt. 0)&
          allocate(rproblemLevel%RmatrixScalar(&
          rproblemDescriptor%nmatrixScalar))

      if (rproblemDescriptor%nmatrixBlock .gt. 0)&
          allocate(rproblemLevel%RmatrixBlock(&
          rproblemDescriptor%nmatrixBlock))

      if (rproblemDescriptor%nvectorScalar .gt. 0)&
          allocate(rproblemLevel%RvectorScalar(&
          rproblemDescriptor%nvectorScalar))

      if (rproblemDescriptor%nvectorBlock .gt. 0)&
          allocate(rproblemLevel%RvectorBlock(&
          rproblemDescriptor%nvectorBlock))
      
#ifdef ENABLE_COPROCESSOR_SUPPORT
      ! Allocate and initialise CUDA streams
      if (rproblemDescriptor%nstream .gt. 0) then
        allocate(rproblemLevel%Istream(&
            rproblemDescriptor%nstream))
        do istream = 1, rproblemDescriptor%nstream
          call coproc_createStream(rproblemLevel%Istream(istream))
        end do
      end if
#endif
      
    end subroutine allocateStructure

  end subroutine problem_initProblemDirect

  !*****************************************************************************

!<subroutine>

  subroutine problem_initProblemFromParlst(rproblem, rparlist, ssectionName)

!<description>
    ! This subroutine initialises the complete problem structure
    ! with the values provided by the parameter list. 
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<output>
    ! problem data structure
    type(t_problem), intent(out) :: rproblem
!</output>
!</subroutine>
    
    ! local variables
    type(t_problemDescriptor) :: rproblemDescriptor

    ! Initialise the problem descriptor
    call problem_initDescriptor(rproblemDescriptor, rparlist, ssectionName)

    ! Initialise the problem rudimentary
    call problem_initProblem(rproblem, rproblemDescriptor)

    ! Now all internal data structures of the gobal problem structure
    ! are prepared. Let us prepare the generation of content and
    ! update the content afterwards. Not that, e.g. matrices depend on
    ! discretisation structures. Hence, the order of the following
    ! calls is crucial to ensure that prerequisite data is generated.
    ! Also note that in the very first update after the general
    ! initialisation, not all internal data may be set correctly so
    ! they need to be updated. Thus, we set bupdateInternals=.true.

    ! 1) Block and spatial Discretisation structures
    call recursiveInitDiscr(rproblem, rproblem, rproblemDescriptor, rparlist)    

    ! 2) Cubature info structures
    call recursiveInitCubInfo(rproblem, rproblem, rproblemDescriptor, rparlist)    

    ! 3) Scalar and block matrices
    call recursiveInitMatrix(rproblem, rproblem, rproblemDescriptor, rparlist)    

    ! 4) Scalar and block vectors
    call recursiveInitVector(rproblem, rproblem, rproblemDescriptor, rparlist)    

    ! 5) Blocks of group finite element structures
    call recursiveInitGroupFEMBlock(rproblem, rproblem, rproblemDescriptor, rparlist)    

    ! 6) Stabilisation structures of AFC-type
    call recursiveInitAFCstab(rproblem, rproblem, rproblemDescriptor, rparlist)    

    ! That`s it
    call problem_releaseDescriptor(rproblemDescriptor)
    
  contains

    ! Here, the real initialisation routine follows.

    !**************************************************************
    ! Internal routine which initialises the discretisation structures
    ! for the entire problem structure
    recursive subroutine recursiveInitDiscr(rproblemTopLevel,&
        rproblem, rproblemDescriptor, rparlist)

      type(t_problem), intent(inout) :: rproblemTopLevel,rproblem
      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
      type(t_parlist), intent(in) :: rparlist
      
      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      integer :: i

      !-------------------------------------------------------------------------
      ! Part 1: Initialise content on all problem levels
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))

        ! Discretisations?
        if (rproblemDescriptor%ndiscretisation .gt. 0) then
          call problem_updateDiscrAll(p_rproblemLevel,&
              rparlist, rproblem%cproblem, rproblemTopLevel%p_rdiscrTasklist,&
              rproblemTopLevel)
        end if
    
        ! Proceed to next problem level
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      !-------------------------------------------------------------------------
      ! Part 2: Initialise content for all subproblems
      if (associated(rproblemDescriptor%Rsubproblem)) then
        p_rproblem => rproblem%p_rproblemFirst
        do i = 1, size(rproblemDescriptor%Rsubproblem)
          
          call recursiveInitDiscr(rproblemTopLevel, p_rproblem,&
              rproblemDescriptor%Rsubproblem(i), rparlist)
          
          ! Proceed to next subproblem
          p_rproblem => p_rproblem%p_rproblemNext
        end do
      end if

    end subroutine recursiveInitDiscr

    !**************************************************************
    ! Internal routine which initialises the cubature info structures
    ! for the entire problem structure
    recursive subroutine recursiveInitCubInfo(rproblemTopLevel,&
        rproblem, rproblemDescriptor, rparlist)

      type(t_problem), intent(inout) :: rproblemTopLevel,rproblem
      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
      type(t_parlist), intent(in) :: rparlist
      
      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      integer :: i

      !-------------------------------------------------------------------------
      ! Part 1: Initialise content on all problem levels
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))

        ! Cubature info structures?
        if (rproblemDescriptor%ncubatureinfo .gt. 0) then
          call problem_updateCubInfoAll(p_rproblemLevel,&
              rparlist, rproblem%cproblem, rproblemTopLevel%p_rcubinfoTasklist,&
              rproblemTopLevel)
        end if

        ! Proceed to next problem level
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      !-------------------------------------------------------------------------
      ! Part 2: Initialise content for all subproblems
      if (associated(rproblemDescriptor%Rsubproblem)) then
        p_rproblem => rproblem%p_rproblemFirst
        do i = 1, size(rproblemDescriptor%Rsubproblem)
          
          call recursiveInitCubInfo(rproblemTopLevel, p_rproblem,&
              rproblemDescriptor%Rsubproblem(i), rparlist)
          
          ! Proceed to next subproblem
          p_rproblem => p_rproblem%p_rproblemNext
        end do
      end if

    end subroutine recursiveInitCubInfo
    
    !**************************************************************
    ! Internal routine which initialises the scalar and block matrices
    ! for the entire problem structure
    recursive subroutine recursiveInitMatrix(rproblemTopLevel,&
        rproblem, rproblemDescriptor, rparlist)

      type(t_problem), intent(inout) :: rproblemTopLevel,rproblem
      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
      type(t_parlist), intent(in) :: rparlist
      
      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      integer :: i

      !-------------------------------------------------------------------------
      ! Part 1: Initialise content on all problem levels
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))
        
        ! Scalar and/or block matrices?
        if ((rproblemDescriptor%nmatrixScalar .gt. 0) .or.&
            (rproblemDescriptor%nmatrixBlock  .gt. 0)) then
          call problem_updateMatrixAll(p_rproblemLevel,&
              rparlist, rproblem%cproblem, rproblemTopLevel%p_rmatrixTasklist,&
              rproblemTopLevel)
        end if
        
        ! Proceed to next problem level
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      !-------------------------------------------------------------------------
      ! Part 2: Initialise content for all subproblems
      if (associated(rproblemDescriptor%Rsubproblem)) then
        p_rproblem => rproblem%p_rproblemFirst
        do i = 1, size(rproblemDescriptor%Rsubproblem)
          
          call recursiveInitMatrix(rproblemTopLevel, p_rproblem,&
              rproblemDescriptor%Rsubproblem(i), rparlist)
          
          ! Proceed to next subproblem
          p_rproblem => p_rproblem%p_rproblemNext
        end do
      end if

    end subroutine recursiveInitMatrix

    !**************************************************************
    ! Internal routine which initialises the scalar and block vectors
    ! for the entire problem structure
    recursive subroutine recursiveInitVector(rproblemTopLevel,&
        rproblem, rproblemDescriptor, rparlist)

      type(t_problem), intent(inout) :: rproblemTopLevel,rproblem
      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
      type(t_parlist), intent(in) :: rparlist
      
      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      integer :: i

      !-------------------------------------------------------------------------
      ! Part 1: Initialise content on all problem levels
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))
        
        ! Scalar and/or block matrices?
        if ((rproblemDescriptor%nvectorScalar .gt. 0) .or.&
            (rproblemDescriptor%nvectorBlock  .gt. 0)) then
          call problem_updateVectorAll(p_rproblemLevel,&
              rparlist, rproblem%cproblem, rproblemTopLevel%p_rvectorTasklist,&
              rproblemTopLevel)
        end if
        
        ! Proceed to next problem level
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      !-------------------------------------------------------------------------
      ! Part 2: Initialise content for all subproblems
      if (associated(rproblemDescriptor%Rsubproblem)) then
        p_rproblem => rproblem%p_rproblemFirst
        do i = 1, size(rproblemDescriptor%Rsubproblem)
          
          call recursiveInitVector(rproblemTopLevel, p_rproblem,&
              rproblemDescriptor%Rsubproblem(i), rparlist)
          
          ! Proceed to next subproblem
          p_rproblem => p_rproblem%p_rproblemNext
        end do
      end if

    end subroutine recursiveInitVector

    !**************************************************************
    ! Internal routine which initialises the blocks of group finite
    ! element structures for the entire problem structure
    recursive subroutine recursiveInitGroupFEMBlock(rproblemTopLevel,&
        rproblem, rproblemDescriptor, rparlist)

      type(t_problem), intent(inout) :: rproblemTopLevel,rproblem
      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
      type(t_parlist), intent(in) :: rparlist
      
      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      integer :: i

      !-------------------------------------------------------------------------
      ! Part 1: Initialise content on all problem levels
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))
        
        ! Scalar and/or block matrices?
        if (rproblemDescriptor%nafcstab .gt. 0) then
          call problem_updateGroupFEMBlockAll(p_rproblemLevel,&
              rparlist, rproblem%cproblem, rproblemTopLevel%p_rgroupFEMBlockTasklist,&
              rproblemTopLevel)
        end if
        
        ! Proceed to next problem level
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      !-------------------------------------------------------------------------
      ! Part 2: Initialise content for all subproblems
      if (associated(rproblemDescriptor%Rsubproblem)) then
        p_rproblem => rproblem%p_rproblemFirst
        do i = 1, size(rproblemDescriptor%Rsubproblem)
          
          call recursiveInitVector(rproblemTopLevel, p_rproblem,&
              rproblemDescriptor%Rsubproblem(i), rparlist)
          
          ! Proceed to next subproblem
          p_rproblem => p_rproblem%p_rproblemNext
        end do
      end if

    end subroutine recursiveInitGroupFEMBlock

    !**************************************************************
    ! Internal routine which initialises the stabilisation structures
    ! of AFC-type for the entire problem structure
    recursive subroutine recursiveInitAFCstab(rproblemTopLevel,&
        rproblem, rproblemDescriptor, rparlist)

      type(t_problem), intent(inout) :: rproblemTopLevel,rproblem
      type(t_problemDescriptor), intent(in) :: rproblemDescriptor
      type(t_parlist), intent(in) :: rparlist
      
      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      integer :: i

      !-------------------------------------------------------------------------
      ! Part 1: Initialise content on all problem levels
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))
        
        ! Scalar and/or block matrices?
        if (rproblemDescriptor%nafcstab .gt. 0) then
          call problem_updateAFCstabAll(p_rproblemLevel,&
              rparlist, rproblem%cproblem, rproblemTopLevel%p_rafcstabTasklist,&
              rproblemTopLevel)
        end if
        
        ! Proceed to next problem level
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      !-------------------------------------------------------------------------
      ! Part 2: Initialise content for all subproblems
      if (associated(rproblemDescriptor%Rsubproblem)) then
        p_rproblem => rproblem%p_rproblemFirst
        do i = 1, size(rproblemDescriptor%Rsubproblem)
          
          call recursiveInitVector(rproblemTopLevel, p_rproblem,&
              rproblemDescriptor%Rsubproblem(i), rparlist)
          
          ! Proceed to next subproblem
          p_rproblem => p_rproblem%p_rproblemNext
        end do
      end if

    end subroutine recursiveInitAFCstab
    
  end subroutine problem_initProblemFromParlst

  !*****************************************************************************

!<subroutine>

  subroutine problem_createProblem(rproblem)

!<description>
    ! This subroutine creates a new problem structure
!</description>

!<intputoutput>
    ! problem data structure
    type(t_problem), intent(out) :: rproblem
!</inputoutput>
!</subroutine>

    ! Reset data
    nullify(rproblem%p_rproblemFirst, rproblem%p_rproblemLast)
    nullify(rproblem%p_rproblemPrev, rproblem%p_rproblemNext)
    nullify(rproblem%p_rproblemLevelMin, rproblem%p_rproblemLevelMax)

  end subroutine problem_createProblem

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_releaseProblem(rproblem, breleaseSubproblems)

!<description>
    ! This subroutine releases an existing problem structure
!</description>

!<input>
    ! OPTIONAL: If present and breleaseSubproblems=TRUE then all
    ! subproblems of the top-level problem rproblem are released.
    logical, optional :: breleaseSubproblems
!</input>

!<inputoutput>
    ! problem data structure
    type(t_problem), intent(inout), target :: rproblem
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_problem), pointer :: p_rproblem
    type(t_problemLevel), pointer :: p_rproblemLevel
    logical :: brelease

    ! Release all subproblems?
    brelease = .false.
    if (present(breleaseSubproblems)) brelease = breleaseSubproblems

    if (brelease) then
      if (associated(rproblem%p_rproblemFirst) .and.&
          associated(rproblem%p_rproblemLast)) then

        ! Start with first subproblem
        p_rproblem => rproblem%p_rproblemFirst

        ! Loop over all subproblems
        do while (associated(p_rproblem))
          call problem_removeProblem(p_rproblem, rproblem)
          call problem_releaseProblem(p_rproblem, breleaseSubproblems)

          deallocate(p_rproblem)
          p_rproblem => rproblem%p_rproblemFirst
        end do
      end if
    end if

    ! Otherwise, we release just the given problem
    p_rproblemLevel => rproblem%p_rproblemLevelMax

    ! Loop over all problem levels and destroy them
    do while(associated(p_rproblemLevel))
      call problem_removeLevel(p_rproblemLevel, rproblem)
      call problem_releaseLevel(p_rproblemLevel)

      deallocate(p_rproblemLevel)
      p_rproblemLevel => rproblem%p_rproblemLevelMax
    end do

    ! Release boundary
    if (associated(rproblem%rboundary)) then
      call boundary_release(rproblem%rboundary)
      deallocate(rproblem%rboundary)
    end if

    ! Release task lists
    if (associated(rproblem%p_rdiscrTasklist))&
        call problem_clearDiscrTL(rproblem%p_rdiscrTasklist)
    if (associated(rproblem%p_rcubinfoTasklist))&
        call problem_clearCubInfoTL(rproblem%p_rcubInfoTasklist)
    if (associated(rproblem%p_rmatrixTasklist))&
        call problem_clearMatrixTL(rproblem%p_rmatrixTasklist)

    ! Reset data
    rproblem%cproblem = ''
    rproblem%cproblemtype = PROBTYPE_UNDEFINED
    nullify(rproblem%p_rproblemTop)
    nullify(rproblem%p_rproblemLevelMin, rproblem%p_rproblemLevelMax)
    nullify(rproblem%p_rproblemPrev, rproblem%p_rproblemNext)
    nullify(rproblem%p_rproblemFirst, rproblem%p_rproblemLast)

  end subroutine problem_releaseProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_appendProblem(rproblem, rsubproblem)

!<description>
    ! This subroutine appends problem structure rsubproblem to the
    ! list of subproblems of problem structure rproblem.
!</description>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! problem structure to be appended
    type(t_problem), intent(inout), target :: rsubproblem
!</inputoutput>
!</subroutine>

    ! Make sure that the subproblem knows its master
    rsubproblem%p_rproblemTop => rproblem

    ! Make sure that the master knows its subproblem
    if (associated(rproblem%p_rproblemFirst) .and.&
        associated(rproblem%p_rproblemLast)) then

      ! Master problem structure has at least one subproblem
      rsubproblem%p_rproblemPrev => rproblem%p_rproblemLast
      nullify(rsubproblem%p_rproblemNext)
      
      rproblem%p_rproblemLast%p_rproblemNext => rsubproblem
      rproblem%p_rproblemLast => rsubproblem

    else 

      ! Problem structure is complete empty
      rproblem%p_rproblemFirst => rsubproblem
      rproblem%p_rproblemLast  => rsubproblem

    end if
    
  end subroutine problem_appendProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_prependProblem(rproblem, rsubproblem)

!<description>
    ! This subroutine prepends problem structure rsubproblem to the
    ! list of subproblems of problem structure rproblem.
!</description>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! problem structure to be prepended
    type(t_problem), intent(inout), target :: rsubproblem
!</inputoutput>
!</subroutine>

    ! Make sure that the subproblem knows its master
    rsubproblem%p_rproblemTop => rproblem

    ! Make sure that the master knows its subproblem
    if (associated(rproblem%p_rproblemFirst) .and.&
        associated(rproblem%p_rproblemLast)) then

      ! Problem structure has at least one subproblem
      rsubproblem%p_rproblemNext => rproblem%p_rproblemFirst
      nullify(rsubproblem%p_rproblemPrev)
      
      rproblem%p_rproblemFirst%p_rproblemPrev => rsubproblem
      rproblem%p_rproblemFirst => rsubproblem

    else 

      ! Problem structure is complete empty
      rproblem%p_rproblemFirst => rsubproblem
      rproblem%p_rproblemLast  => rsubproblem

    end if

  end subroutine problem_prependProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_removeProblem(rproblem, rproblemRef)

!<description>
    ! This subroutine removes an existing problem structure from the
    ! linked list of problem structures. If the optional reference
    ! problem structure rproblemRef is provided it is assumed that
    ! rproblem is a subproblem of rproblemRef.
!</description>

!<inputoutput>
    ! problem structure to be removed
    type(t_problem), intent(inout) :: rproblem

    ! OPTIONAL: reference problem structure
    type(t_problem), intent(inout), target, optional :: rproblemRef
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_problem), pointer :: p_rproblem

    ! Do we have a reference problem?
    p_rproblem => rproblem%p_rproblemTop
    if (present(rproblemRef)) p_rproblem => rproblemRef

    ! Do we have a top-level problem? If not, then we are done.
    if (associated(rproblem%p_rproblemTop)) then
      
      ! Check consistency
      if (.not.associated(rproblem%p_rproblemTop, p_rproblem)) then
        call output_line('Problem structure does not belong to top-level problem structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_removeProblem')
        call sys_halt()
      end if

      ! Is this the first problem?
      if (.not.associated(rproblem%p_rproblemPrev)) then
        p_rproblem%p_rproblemFirst => rproblem%p_rproblemNext
      else
        rproblem%p_rproblemPrev%p_rproblemNext => rproblem%p_rproblemNext
      end if
      
      ! Is this the last problem?
      if (.not.associated(rproblem%p_rproblemNext)) then
        p_rproblem%p_rproblemLast => rproblem%p_rproblemPrev
      else
        rproblem%p_rproblemNext%p_rproblemPrev => rproblem%p_rproblemPrev
      end if
    end if
      
  end subroutine problem_removeProblem

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_infoProblem(rproblem, bverbose)

!<description>
    ! This subroutine outputs information about the problem structure
!</description>

!<input>
    ! problem structure
    type(t_problem), intent(in) :: rproblem

    ! OPTIONAL: switch to turn on verbose output
    ! If not present or .FALSE. then standard output is used
    logical, intent(in), optional :: bverbose
!</input>
!</subroutine>

    ! local variable
    logical :: bverb
    
    bverb = .false.
    if (present(bverbose)) bverb = bverbose

    if (bverb) then
      ! Verbose output
      call doInfoVerbose(rproblem)
    else
      ! Standard output
      call doInfoStandard(rproblem, 0)
    end if

  contains

    ! Here, the real info routines follow.

    !**************************************************************
    ! Internal routine which outputs verbose information about the
    ! entire problem structure including all subproblems
    recursive subroutine doInfoVerbose(rproblem)
      
      type(t_problem), intent(in) :: rproblem

      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      integer :: i

      call output_lbrk()
      call output_separator(OU_SEP_HASH)
      call output_line('Problem: '//trim(rproblem%cproblem))
      call output_separator(OU_SEP_HASH)
      call output_lbrk()
      
      !-------------------------------------------------------------------------
      ! Part 1: Display content on all problem levels
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))

        call output_separator(OU_SEP_MINUS)
        call output_line('Problemlevel: '//trim(sys_siL(p_rproblemLevel%ilev,3)))
        
        ! Discretisations
        if (associated(p_rproblemLevel%Rdiscretisation)) then
          do i=1,size(p_rproblemLevel%Rdiscretisation)
            call output_separator(OU_SEP_TILDE)
            call output_line('Blockdiscretisation: '//trim(sys_siL(i,3)))
            call spdiscr_infoBlockDiscr(p_rproblemLevel%Rdiscretisation(i))
          end do
        end if

        ! Cubature info structures
        if (associated(p_rproblemLevel%RcubatureInfo)) then
          do i=1,size(p_rproblemLevel%RcubatureInfo)
            call output_separator(OU_SEP_TILDE)
            call output_line('Cubature info structure: '//trim(sys_siL(i,3)))
            call spdiscr_infoCubatureInfo(p_rproblemLevel%RcubatureInfo(i))
          end do
        end if

        ! Scalar matrices
        if (associated(p_rproblemLevel%RmatrixScalar)) then
          do i=1,size(p_rproblemLevel%RmatrixScalar)
            call output_separator(OU_SEP_TILDE)
            call output_line('Scalar matrix: '//trim(sys_siL(i,3)))
            call lsyssc_infoMatrix(p_rproblemLevel%RmatrixScalar(i))
          end do
        end if

        ! Block matrices
        if (associated(p_rproblemLevel%RmatrixBlock)) then
          do i=1,size(p_rproblemLevel%RmatrixBlock)
            call output_separator(OU_SEP_TILDE)
            call output_line('Block matrix: '//trim(sys_siL(i,3)))
            call lsysbl_infoMatrix(p_rproblemLevel%RmatrixBlock(i))
          end do
        end if

        ! Scalar vectors
        if (associated(p_rproblemLevel%RvectorScalar)) then
          do i=1,size(p_rproblemLevel%RvectorScalar)
            call output_separator(OU_SEP_TILDE)
            call output_line('Scalar vector: '//trim(sys_siL(i,3)))
            call lsyssc_infoVector(p_rproblemLevel%RvectorScalar(i))
          end do
        end if

        ! Block vectors
        if (associated(p_rproblemLevel%RvectorBlock)) then
          do i=1,size(p_rproblemLevel%RvectorBlock)
            call output_separator(OU_SEP_TILDE)
            call output_line('Block vector: '//trim(sys_siL(i,3)))
            call lsysbl_infoVector(p_rproblemLevel%RvectorBlock(i))
          end do
        end if

        ! Block of group finite element structures
        if (associated(p_rproblemLevel%RgroupFEMBlock)) then
          do i=1,size(p_rproblemLevel%RgroupFEMBlock)
            call output_separator(OU_SEP_TILDE)
            call output_line('Block of group finite element structures: '//trim(sys_siL(i,3)))
            call gfem_infoGroupFEMBlock(p_rproblemLevel%RgroupFEMBlock(i))
          end do
        end if

        ! Stabilisation structures of AFC-type
        if (associated(p_rproblemLevel%Rafcstab)) then
          do i=1,size(p_rproblemLevel%Rafcstab)
            call output_separator(OU_SEP_TILDE)
            call output_line('Stabilisation structure: '//trim(sys_siL(i,3)))
            call afcstab_infoStabilisation(p_rproblemLevel%Rafcstab(i))
          end do
        end if

        ! Proceed to next problem level
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      !-------------------------------------------------------------------------
      ! Part 2: Display content for all subproblems
      p_rproblem => rproblem%p_rproblemFirst
      do while(associated(p_rproblem))
        call doInfoVerbose(p_rproblem)
        
        ! Proceed to next subproblem
        p_rproblem => p_rproblem%p_rproblemNext
      end do
        
    end subroutine doInfoVerbose

    !**************************************************************

    recursive subroutine doInfoStandard(rproblem, nindent)

      type(t_problem), intent(in), target :: rproblem
      integer, intent(in) :: nindent

      ! local variables
      type(t_problem), pointer :: p_rproblem
      type(t_problemLevel), pointer :: p_rproblemLevel
      character(len=SYS_STRLEN) :: string
      character(len=nindent) :: cindent
      integer :: i

      ! Prepare indentation
      do i=1,nindent
        cindent(i:i) = ' '
      end do
      
      ! Output header
      call output_line (cindent//'Problem            : '//trim(rproblem%cproblem))
      call output_line (cindent//'--------------------')
      
      ! Output problem type(s)
      string = ''
      if (iand(rproblem%cproblemtype, PROBTYPE_META) .ne. 0)&
          string = trim(string)//' META'
      if (iand(rproblem%cproblemtype, PROBTYPE_STANDALONE) .ne. 0)&
          string = trim(string)//' STANDALONE'
      if (iand(rproblem%cproblemtype, PROBTYPE_EMBEDDED) .ne. 0)&
          string = trim(string)//' EMBEDDED'
      call output_line (cindent//'Problem type       :'//trim(string))

      if (associated(rproblem%p_rproblemTop)) then
        call output_line (cindent//'Top-level problem  : '//&
            trim(rproblem%p_rproblemTop%cproblem))
      end if

      if (associated(rproblem%p_rproblemFirst)) then
        call output_line (cindent//'List of subproblems')
        call output_line (cindent//'-------------------')

        p_rproblem => rproblem%p_rproblemFirst
        string = ''
        do while (associated(p_rproblem))
          string = trim(string)//' '//trim(p_rproblem%cproblem)
          if (associated(p_rproblem, rproblem%p_rproblemLast)) exit
          p_rproblem => p_rproblem%p_rproblemNext
        end do
        call output_line (cindent//'Subproblem         :'//trim(string))

      end if

      if (associated(rproblem%p_rproblemLevelMin)) then
        call output_line (cindent//'minimum level      : '//&
            trim(sys_siL(rproblem%p_rproblemLevelMin%ilev,15)))
      else
        call output_line (cindent//'minimum level      : not associated')
      end if

      if (associated(rproblem%p_rproblemLevelMax)) then
        call output_line (cindent//'maximum level      : '//&
            trim(sys_siL(rproblem%p_rproblemLevelMax%ilev,15)))
      else
        call output_line (cindent//'maximum level      : not associated')
      end if

      call output_lbrk()

      ! Initialisation
      p_rproblemLevel => rproblem%p_rproblemLevelMax

      ! Loop over all problem levels
      do while(associated(p_rproblemLevel))
        call problem_infoLevel(p_rproblemLevel, nindent)
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

      ! Output information about possible subproblems
      if (associated(rproblem%p_rproblemFirst) .and.&
          associated(rproblem%p_rproblemLast)) then
        
        p_rproblem => rproblem%p_rproblemFirst
        do while (associated(p_rproblem))
          call doInfoStandard(p_rproblem, nindent+2)
          if (associated(p_rproblem, rproblem%p_rproblemLast)) exit
          p_rproblem => p_rproblem%p_rproblemNext
        end do
      end if
      
    end subroutine doInfoStandard

  end subroutine problem_infoProblem

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_infoProblemDescriptor(rproblemDescriptor, nindent)

!<description>
    ! This subroutine outputs information about the problem descriptor
!</description>

!<input>
    ! problem descriptor
    type(t_problemDescriptor), intent(in) :: rproblemDescriptor

    ! OPTIONAL: number of indentation white spaces
    integer, intent(in), optional :: nindent
!</input>
!</subroutine>

    ! local variable
    integer :: i,j
    
    ! Prepare indentation
    i = 0
    if (present(nindent)) i=nindent
    
    ! Output information about the current problem
    call doInfo(rproblemDescriptor, i)
    
    ! Output information about possible subproblems
    if (associated(rproblemDescriptor%Rsubproblem)) then
      do j=lbound(rproblemDescriptor%Rsubproblem,1),&
           ubound(rproblemDescriptor%Rsubproblem,1)
        call problem_infoProblemDescriptor(&
            rproblemDescriptor%Rsubproblem(j), i+2)
      end do
    end if
    
  contains

    ! Here, the real info routines follow.
    
    !**************************************************************
    
    subroutine doInfo(rproblemDescriptor, nindent)

      type(t_problemDescriptor), intent(in), target :: rproblemDescriptor
      integer, intent(in) :: nindent

      ! local variables
      character(len=SYS_STRLEN) :: string
      character(len=nindent) :: cindent
      integer :: i

      ! Prepare indentation
      do i=1,nindent
        cindent(i:i) = ' '
      end do
      
      ! Output header
      call output_line (cindent//'Problem            : '//trim(rproblemDescriptor%cproblem))
      call output_line (cindent//'--------------------')

      ! Output problem type(s)
      string = ''
      if (iand(rproblemDescriptor%cproblemtype, PROBTYPE_META) .ne. 0)&
          string = trim(string)//' META'
      if (iand(rproblemDescriptor%cproblemtype, PROBTYPE_STANDALONE) .ne. 0)&
          string = trim(string)//' STANDALONE'
      if (iand(rproblemDescriptor%cproblemtype, PROBTYPE_EMBEDDED) .ne. 0)&
          string = trim(string)//' EMBEDDED'
      call output_line (cindent//'Problem type       :'//trim(string))

      ! Output problem description
      call output_line (cindent//'dimension          : '//&
          trim(sys_siL(rproblemDescriptor%ndimension,15)))
      call output_line (cindent//'minimum level      : '//&
          trim(sys_siL(rproblemDescriptor%nlmin,15)))
      call output_line (cindent//'maximum level      : '//&
          trim(sys_siL(rproblemDescriptor%nlmax,15)))
      call output_line (cindent//'ndiscretisation    : '//&
          trim(sys_siL(rproblemDescriptor%ndiscretisation,15)))
      call output_line (cindent//'ncubatureinfo      : '//&
          trim(sys_siL(rproblemDescriptor%ncubatureInfo,15)))
      call output_line (cindent//'nafcstab           : '//&
          trim(sys_siL(rproblemDescriptor%nafcstab,15)))
      call output_line (cindent//'ngroupfemblock     : '//&
          trim(sys_siL(rproblemDescriptor%ngroupfemBlock,15)))
      call output_line (cindent//'nmatrixscalar      : '//&
          trim(sys_siL(rproblemDescriptor%nmatrixScalar,15)))
      call output_line (cindent//'nmatrixblock       : '//&
          trim(sys_siL(rproblemDescriptor%nmatrixBlock,15)))
      call output_line (cindent//'nvectorscalar      : '//&
          trim(sys_siL(rproblemDescriptor%nvectorScalar,15)))
      call output_line (cindent//'nvectorblock       : '//&
          trim(sys_siL(rproblemDescriptor%nvectorBlock,15)))
#ifdef ENABLE_COPROCESSOR_SUPPORT
      call output_line (cindent//'nstream            : '//&
          trim(sys_siL(rproblemDescriptor%nstream,15)))
#endif
      call output_line (cindent//'trifile            : '//&
          trim(rproblemDescriptor%trifile))
      call output_line (cindent//'prmfile            : '//&
          trim(rproblemDescriptor%prmfile))
      call output_line (cindent//'iconvStrategy      : '//&
          trim(sys_siL(rproblemDescriptor%iconvStrategy,15)))
      call output_line (cindent//'ddisturbMesh       : '//&
          trim(sys_sdL(rproblemDescriptor%ddisturbmesh,15)))
  
    end subroutine doInfo
    
  end subroutine problem_infoProblemDescriptor

  !*****************************************************************************

!<function>

  function problem_getProblem(rproblem, cproblem) result(p_rproblem)

!<description>
    ! This function returns a pointer to the problem structure which
    ! is labeled cproblem. If such problem structure does not exist
    ! p_rproblem points to null().
!</description>

!<input>
    ! top-level problem structure
    type(t_problem), intent(in), target :: rproblem

    ! label of the problem structure to be returned
    character(len=*), intent(in) :: cproblem
!</input>

!<result>
    ! problem structure
    type(t_problem), pointer :: p_rproblem
!</result>
!</function>

    p_rproblem => getProblem(rproblem, cproblem)
    if (associated(p_rproblem)) return

    ! If we are here, then we need to repeat the search starting at
    ! the top most problem structure
    p_rproblem => rproblem
    do while(associated(p_rproblem%p_rproblemTop))
      p_rproblem => p_rproblem%p_rproblemTop
    end do

    ! Ok, so let us repeat the search...
    p_rproblem => getProblem(p_rproblem, cproblem)

  contains

    ! Here, the real working routines follow
    
    !**************************************************************

    recursive function getProblem(rproblem, cproblem) result(p_rproblem)

      type(t_problem), intent(in), target :: rproblem
      character(len=*), intent(in) :: cproblem
      type(t_problem), pointer :: p_rproblem

      ! local variable
      type(t_problem), pointer :: p_rsubproblem
      
      if (trim(rproblem%cproblem) .eq. cproblem) then
        ! That`s it
        p_rproblem => rproblem
        return
      else
        nullify(p_rproblem)
        
        if (associated(rproblem%p_rproblemFirst) .and.&
            associated(rproblem%p_rproblemLast)) then
          ! Let us proceed in all subproblems if they exists        
          p_rsubproblem => rproblem%p_rproblemFirst
          do while (associated(p_rsubproblem))
            p_rproblem => getProblem(p_rsubproblem, cproblem)
            if (associated(p_rproblem)) return
            p_rsubproblem => p_rsubproblem%p_rproblemNext
          end do
        end if
      end if
    end function getProblem

  end function problem_getProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_createLevel(rproblemLevel, ilev)

!<description>
    ! This subroutine creates a new problem level structure
!</description>

!<input>
    ! level number
    integer, intent(in) :: ilev
!</input>

!<output>
    ! problem level structure
    type(t_problemLevel), intent(out) :: rproblemLevel
!</output>
!</subroutine>

    ! Set problem level
    rproblemLevel%ilev = ilev

  end subroutine problem_createLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_releaseLevel(rproblemLevel)

!<description>
    ! This subroutine releases an existing problem level structure
!</description>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i
    
    ! Release triangulation structure?
    if ((iand(rproblemLevel%p_rproblem%cproblemtype,PROBTYPE_META) .eq. 0).and.&
        (iand(rproblemLevel%p_rproblem%cproblemtype,PROBTYPE_EMBEDDED) .eq. 0)) then
      call tria_done(rproblemLevel%rtriangulation)
      deallocate(rproblemLevel%rtriangulation)
    else
      nullify(rproblemLevel%rtriangulation)
    end if

    ! Release discretisation structures
    if (associated(rproblemLevel%Rdiscretisation)) then
      do i = lbound(rproblemLevel%Rdiscretisation,1),&
             ubound(rproblemLevel%Rdiscretisation,1)
        call spdiscr_releaseBlockDiscr(rproblemLevel%Rdiscretisation(i))
      end do
      deallocate(rproblemLevel%Rdiscretisation)
    end if

    ! Release scalar cubature info structures
    if (associated(rproblemLevel%RcubatureInfo)) then
      do i = lbound(rproblemLevel%RcubatureInfo,1),&
             ubound(rproblemLevel%RcubatureInfo,1)
        call spdiscr_releaseCubStructure(rproblemLevel%RcubatureInfo(i))
      end do
      deallocate(rproblemLevel%RcubatureInfo)
    end if

    ! Release all scalar matrices
    if (associated(rproblemLevel%RmatrixScalar)) then
      do i = lbound(rproblemLevel%RmatrixScalar,1),&
             ubound(rproblemLevel%RmatrixScalar,1)
        call lsyssc_releaseMatrix(rproblemLevel%RmatrixScalar(i))
      end do
      deallocate(rproblemLevel%RmatrixScalar)
    end if

    ! Release all block matrices
    if (associated(rproblemLevel%RmatrixBlock)) then
      do i = lbound(rproblemLevel%RmatrixBlock,1),&
             ubound(rproblemLevel%RmatrixBlock,1)
        call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(i))
      end do
      deallocate(rproblemLevel%RmatrixBlock)
    end if

    ! Release all scalar vectors
    if (associated(rproblemLevel%RvectorScalar)) then
      do i = lbound(rproblemLevel%RvectorScalar,1),&
             ubound(rproblemLevel%RvectorScalar,1)
        call lsyssc_releaseVector(rproblemLevel%RvectorScalar(i))
      end do
      deallocate(rproblemLevel%RvectorScalar)
    end if

    ! Release all block vectors
    if (associated(rproblemLevel%RvectorBlock)) then
      do i = lbound(rproblemLevel%RvectorBlock,1),&
             ubound(rproblemLevel%RvectorBlock,1)
        call lsysbl_releaseVector(rproblemLevel%RvectorBlock(i))
      end do
      deallocate(rproblemLevel%RvectorBlock)
    end if

    ! Release stabilisation structure
    if (associated(rproblemLevel%Rafcstab)) then
      do i = lbound(rproblemLevel%Rafcstab,1),&
             ubound(rproblemLevel%Rafcstab,1)
        call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
      end do
      deallocate(rproblemLevel%Rafcstab)
    end if

    ! Release group finite element block
    if (associated(rproblemLevel%RgroupFEMBlock)) then
      do i = lbound(rproblemLevel%RgroupFEMBlock,1),&
             ubound(rproblemLevel%RgroupFEMBlock,1)
        call gfem_releaseGroupFEMBlock(rproblemLevel%RgroupFEMBlock(i))
      end do
      deallocate(rproblemLevel%RgroupFEMBlock)
    end if

#ifdef ENABLE_COPROCESSOR_SUPPORT
    ! Release CUDA streams
    if (associated(rproblemLevel%Istream)) then
      do i = lbound(rproblemLevel%Istream,1),&
             ubound(rproblemLevel%Istream,1)
        call coproc_synchronizeStream(rproblemLevel%Istream(i))
        call coproc_destroyStream(rproblemLevel%Istream(i))
      end do
      deallocate(rproblemLevel%Istream)
    end if
#endif

  end subroutine problem_releaseLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_appendLevel(rproblem, rproblemLevel, rproblemLevelRef)

!<description>
    ! This subroutine appends a problem level structure to the linked
    ! list of problem level structures. If the optional reference
    ! level rproblemLevelRef is given, the problem level structure is
    ! appended to rproblemLevelRef. Otherwise, the maximum problem
    ! level is used as reference structure.
!</description>

!<inputoutput>
    ! global problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! OPTIONAL: reference problem structure
    type(t_problemLevel), intent(inout), target, optional :: rproblemLevelRef
!</inputoutput>
!</subroutine>

    ! Set pointer to global problem structure
    rproblemLevel%p_rproblem => rproblem

    if (associated(rproblem%p_rproblemLevelMin) .and.&
        associated(rproblem%p_rproblemLevelMax)) then

      if (present(rproblemLevelRef)) then

        ! Insert rproblemLevel after rproblemLevelRef
        rproblemLevel%p_rproblemLevelCoarse => rproblemLevelRef
        rproblemLevel%p_rproblemLevelFine   => rproblemLevelRef%p_rproblemLevelFine

        if (associated(rproblemLevelRef%p_rproblemLevelFine)) then
          rproblemLevelRef%p_rproblemLevelFine%p_rproblemLevelCoarse => rproblemLevel
        else
          rproblem%p_rproblemLevelMax => rproblemLevel
        end if

        rproblemLevelRef%p_rproblemLevelFine => rproblemLevel

      else

        ! Set pointer to maximum problem level
        rproblemLevel%p_rproblemLevelCoarse => rproblem%p_rproblemLevelMax
        nullify(rproblemLevel%p_rproblemLevelFine)

        rproblem%p_rproblemLevelMax%p_rproblemLevelFine => rproblemLevel
        rproblem%p_rproblemLevelMax => rproblemLevel

      end if

    else

      ! Problem structure is completely empty
      rproblem%p_rproblemLevelMin => rproblemLevel
      rproblem%p_rproblemLevelMax => rproblemLevel

      nullify(rproblemLevel%p_rproblemLevelCoarse, rproblemLevel%p_rproblemLevelFine)

    end if
  end subroutine problem_appendLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_prependLevel(rproblem, rproblemLevel, rproblemLevelRef)

!<description>
    ! This subroutine prepends a problem level structure to the linked
    ! list of problem level structures. If the optional reference
    ! level rproblemLevelRef is given, the problem level structure is
    ! prepended to rproblemLevelRef. Otherwise, the minimum problem
    ! level is used as reference structure.
!</description>

!<inputoutput>
    ! global problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! OPTIONAL: reference problem structure
    type(t_problemLevel), intent(inout), target, optional :: rproblemLevelRef
!</inputoutput>
!</subroutine>

    ! Set pointer to global problem structure
    rproblemLevel%p_rproblem => rproblem

    if (associated(rproblem%p_rproblemLevelMin) .and.&
        associated(rproblem%p_rproblemLevelMax)) then

      if (present(rproblemLevelRef)) then

        ! Insert rproblemLevel before rproblemLevelRef
        rproblemLevel%p_rproblemLevelCoarse => rproblemLevelRef%p_rproblemLevelCoarse
        rproblemLevel%p_rproblemLevelFine   => rproblemLevelRef

        if (associated(rproblemLevelRef%p_rproblemLevelCoarse)) then
          rproblemLevelRef%p_rproblemLevelCoarse%p_rproblemLevelFine => rproblemLevel
        else
          rproblem%p_rproblemLevelMin => rproblemLevel
        end if

        rproblemLevelRef%p_rproblemLevelCoarse => rproblemLevel

      else

        ! Set pointer to minimum problem level
        rproblemLevel%p_rproblemLevelFine => rproblem%p_rproblemLevelMin
        nullify(rproblemLevel%p_rproblemLevelCoarse)

        rproblem%p_rproblemLevelMin%p_rproblemLevelCoarse => rproblemLevel
        rproblem%p_rproblemLevelMin => rproblemLevel

      end if

    else

      ! Problem structure is completely empty
      rproblem%p_rproblemLevelMin => rproblemLevel
      rproblem%p_rproblemLevelMax => rproblemLevel

      nullify(rproblemLevel%p_rproblemLevelCoarse, rproblemLevel%p_rproblemLevelFine)

    end if

  end subroutine problem_prependLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_removeLevel(rproblemLevel, rproblemRef)

!<description>
    ! This subroutine removes an existing problem level structure from
    ! an existing problem structure
!</description>

!<inputoutput>
    ! problem level structure to be removed
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! OPTIONAL: reference problem structure
    type(t_problem), intent(inout), target, optional :: rproblemRef
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_problem), pointer :: p_rproblem

    ! Do we have a reference problem?
    p_rproblem => rproblemLevel%p_rproblem
    if (present(rproblemRef)) p_rproblem => rproblemRef

    ! Check consistency
    if (.not.associated(rproblemLevel%p_rproblem, p_rproblem)) then
      call output_line('Problem level structure does not belong to problem structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_removeLevel')
      call sys_halt()
    end if

    ! Is this the coarsest problem?
    if (.not.associated(rproblemLevel%p_rproblemLevelCoarse)) then
      p_rproblem%p_rproblemLevelMin => rproblemLevel%p_rproblemLevelFine
    else
      rproblemLevel%p_rproblemLevelCoarse%p_rproblemLevelFine => rproblemLevel%p_rproblemLevelFine
    end if

    ! Is this the finest problem?
    if (.not.associated(rproblemLevel%p_rproblemLevelFine)) then
      p_rproblem%p_rproblemLevelMax => rproblemLevel%p_rproblemLevelCoarse
    else
      rproblemLevel%p_rproblemLevelFine%p_rproblemLevelCoarse => rproblemLevel%p_rproblemLevelCoarse
    end if

  end subroutine problem_removeLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_infoLevel(rproblemLevel, nindent)

!<description>
    ! This subroutine outputs information about the problem level structure
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! OPTIONAL: number of indentation white spaces
    integer, intent(in), optional :: nindent
!</input>
!</subroutine>

    ! local variable
    integer :: i

    ! Prepare indentation
    i = 0
    if (present(nindent)) i=nindent
    
    call doInfo(rproblemLevel, i)

  contains
    
    ! Here, the real info routines follow.
    
    !**************************************************************
    
    subroutine doInfo(rproblemLevel, nindent)

      type(t_problemLevel), intent(in), target :: rproblemLevel
      integer, intent(in) :: nindent

      ! local variables
      character(len=nindent) :: cindent
      integer :: i
      
      ! Prepare indentation
      do i=1,nindent
        cindent(i:i) = ' '
      end do
      
      ! Output information about the current level
      call output_line (cindent//'ProblemLevel       : '//trim(sys_siL(rproblemLevel%ilev,15)))
      call output_line (cindent//'--------------------')
      
      if (associated(rproblemLevel%p_rproblemLevelCoarse)) then
        call output_line (cindent//'coarser level      : '//&
                          trim(sys_siL(rproblemLevel%p_rproblemLevelCoarse%ilev,15)))
      else
        call output_line (cindent//'coarser level      : not associated')
      end if
      
      if (associated(rproblemLevel%p_rproblemLevelFine)) then
        call output_line (cindent//'finer level        : '//&
                          trim(sys_siL(rproblemLevel%p_rproblemLevelFine%ilev,15)))
      else
        call output_line (cindent//'finer level        : not associated')
      end if
      
      if (associated(rproblemLevel%Rafcstab)) then
        call output_line (cindent//'Rafcstab           : '//&
                          trim(sys_siL(size(rproblemLevel%Rafcstab),15)))
      else
        call output_line (cindent//'Rafcstab           : not associated')
      end if
      
      if (associated(rproblemLevel%RgroupFEMBlock)) then
        call output_line (cindent//'RgroupFEMBlock     : '//&
                          trim(sys_siL(size(rproblemLevel%RgroupFEMBlock),15)))
      else
        call output_line (cindent//'RgroupFEMBlock     : not associated')
      end if

      if (associated(rproblemLevel%RmatrixScalar)) then
        call output_line (cindent//'RmatrixScalar      : '//&
                          trim(sys_siL(size(rproblemLevel%RmatrixScalar),15)))
      else
        call output_line (cindent//'RmatrixScalar      : not associated')
      end if
      
      if (associated(rproblemLevel%RmatrixBlock)) then
        call output_line (cindent//'RmatrixBlock       : '//&
                          trim(sys_siL(size(rproblemLevel%RmatrixBlock),15)))
      else
        call output_line (cindent//'RmatrixBlock       : not associated')
      end if

      if (associated(rproblemLevel%RvectorScalar)) then
        call output_line (cindent//'RvectorScalar      : '//&
                         trim(sys_siL(size(rproblemLevel%RvectorScalar),15)))
      else
        call output_line (cindent//'RvectorScalar      : not associated')
      end if

      if (associated(rproblemLevel%RvectorBlock)) then
        call output_line (cindent//'RvectorBlock       : '//&
                          trim(sys_siL(size(rproblemLevel%RvectorBlock),15)))
      else
        call output_line (cindent//'RvectorBlock       : not associated')
      end if
      call output_lbrk()

    end subroutine doInfo

  end subroutine problem_infoLevel

  !*****************************************************************************

!<function>

  function problem_getLevelDirect(rproblem, ilev, btopdown)&
      result(p_rproblemLevel)

!<description>
    ! This subroutine returns a pointer to the problem level structure
    ! with specified level number ilev. If such problem level does not
    ! exist p_rproblemLevel points to null().
!</description>

!<input>
    ! global problem structure
    type(t_problem), intent(in) :: rproblem

    ! level number
    integer, intent(in) :: ilev

    ! OPTIONAL: search direction
    ! If this parameter is .false., the loop starts at the minimum
    ! problem level and proceeds to the next finer problem
    ! level. Otherwise, the loop starts at the maximum problem level
    ! and continues with the next coarser level. The default value is
    ! .true., that is, top-down search is performed
    logical, intent(in), optional :: btopdown
!</input>

!<result
    ! problem level structure
    type(t_problemLevel), pointer :: p_rproblemLevel
!</result>
!</function>

    ! local variable
    logical :: bisTopdown

    bisTopdown = .true.
    if (present(btopdown)) bisTopdown=btopdown

    if (bisTopdown) then

      ! Top-down search
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while (associated(p_rproblemLevel))
        if(p_rproblemLevel%ilev .eq. ilev) return
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

    else

      ! Bottom-up search
      p_rproblemLevel => rproblem%p_rproblemLevelMin
      do while (associated(p_rproblemLevel))
        if(p_rproblemLevel%ilev .eq. ilev) return
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelFine
      end do

    end if

    nullify(p_rproblemLevel)

  end function problem_getLevelDirect

  !*****************************************************************************

!<function>

  function problem_getLevelIndirect(rproblem, cproblem, ilev, btopdown)&
      result(p_rproblemLevel)

!<description>
    ! This subroutine returns a pointer to the problem level structure
    ! with specified level number ilev in the problem structure
    ! labeled cproblem. If such problem (level) does not exist
    ! p_rproblemLevel points to null().
!</description>

!<input>
    ! top-level problem structure
    type(t_problem), intent(in), target :: rproblem
    
    ! label of the problem structure to be returned
    character(len=*), intent(in) :: cproblem

    ! level number
    integer, intent(in) :: ilev

    ! OPTIONAL: search direction
    ! If this parameter is .false., the loop starts at the minimum
    ! problem level and proceeds to the next finer problem
    ! level. Otherwise, the loop starts at the maximum problem level
    ! and continues with the next coarser level. The default value is
    ! .true., that is, top-down search is performed
    logical, intent(in), optional :: btopdown
!</input>

!<result
    ! problem level structure
    type(t_problemLevel), pointer :: p_rproblemLevel
!</result>
!</function>

    ! local variable
    type(t_problem), pointer :: p_rproblem

    ! Find problem structure by label
    p_rproblem => problem_getProblem(rproblem, cproblem)
    if (associated(p_rproblem)) then
      p_rproblemLevel => problem_getLevel(p_rproblem, ilev, btopdown)
    else
      nullify(p_rproblemLevel)
    end if

  end function problem_getLevelIndirect

  !*****************************************************************************

!<subroutine>

  subroutine problem_setSpec(rproblem, iSpec, coperation, ilev)

!<description>
    ! This subroutine combines the problem level specification flag
    ! with the bitfield iSpec using one of the operations (iand, ior,
    ! ieor) specified by ioperation.
!</description>

!<input>
    ! Bit field to set
    integer(I32), intent(in) :: iSpec

    ! Operation to combine iSpec and the problem level specification flag
    character(len=*), intent(in) :: coperation

    ! OPTIONAL: level number
    integer, intent(in), optional :: ilev
!</input>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! local variable
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    
    if (present(ilev)) then
      
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))
        if (p_rproblemLevel%ilev .eq. ilev) then
          if ((coperation .eq. 'iand') .or.&
              (coperation .eq. 'IAND')) then
            p_rproblemLevel%iproblemSpec = iand(p_rproblemLevel%iproblemSpec, iSpec)
          elseif ((coperation .eq. 'ior') .or.&
                  (coperation .eq. 'IOR')) then
            p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec, iSpec)
          elseif ((coperation .eq. 'ieor') .or.&
                  (coperation .eq. 'IEOR')) then
            p_rproblemLevel%iproblemSpec = ieor(p_rproblemLevel%iproblemSpec, iSpec)
          else
            p_rproblemLevel%iproblemSpec = iSpec
          end if
          return
        end if
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

    else
      
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while(associated(p_rproblemLevel))
        if ((coperation .eq. 'iand') .or.&
            (coperation .eq. 'IAND')) then
          p_rproblemLevel%iproblemSpec = iand(p_rproblemLevel%iproblemSpec, iSpec)
        elseif ((coperation .eq. 'ior') .or.&
                (coperation .eq. 'IOR')) then
          p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec, iSpec)
        elseif ((coperation .eq. 'ieor') .or.&
                (coperation .eq. 'IEOR')) then
          p_rproblemLevel%iproblemSpec = ieor(p_rproblemLevel%iproblemSpec, iSpec)
        else
          p_rproblemLevel%iproblemSpec = iSpec
        end if
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

    end if

  end subroutine problem_setSpec

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_initDescriptor(rproblemDescriptor, rparlist, ssectionName)

!<description>
    ! This subroutine initialises the problem descriptor of a complete
    ! problem structure with the values provided by the parameter list.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<output>
    ! problem descriptor
    type(t_problemDescriptor), intent(out) :: rproblemDescriptor
!</output>
!</subroutine>
    
    ! local variables
    character(len=SYS_STRLEN) :: sproblemType,sproblemName
    integer :: i,nsubproblem

    ! Let us create the problem descriptor
    rproblemDescriptor%cproblem = ssectionName
    call sys_toupper(rproblemDescriptor%cproblem)

    ! Get problem type
    call parlst_getvalue_string(rparlist, ssectionName,&
        'problemtype', sproblemType, bdequote=.true.)
    call sys_toupper(sproblemType)

    ! What type of problem are we?
    if (trim(sproblemType) .eq. 'META') then
      ! This is a meta-problem
      rproblemDescriptor%cproblemtype = ior(rproblemDescriptor%cproblemtype,&
                                            PROBTYPE_META)
    elseif (trim(sproblemType) .eq. 'STANDALONE') then
      ! This is a standalone problem
      rproblemDescriptor%cproblemtype = ior(rproblemDescriptor%cproblemtype,&
                                            PROBTYPE_STANDALONE)
    elseif (trim(sproblemType) .eq. 'EMBEDDED') then
      ! This is an embedded problem
      rproblemDescriptor%cproblemtype = ior(rproblemDescriptor%cproblemtype,&
                                            PROBTYPE_EMBEDDED)
      ! Moreover, it is also a standalone problem
      rproblemDescriptor%cproblemtype = ior(rproblemDescriptor%cproblemtype,&
                                            PROBTYPE_STANDALONE)
    else
      call output_line('Unsupported problem type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_initDescriptor')
      call sys_halt()
    end if

    ! Do we have one or more subproblems?
    nsubproblem = parlst_querysubstrings(rparlist, ssectionName, 'subproblem')
    if (nsubproblem .gt. 0) then
      allocate(rproblemDescriptor%Rsubproblem(nsubproblem))
      
      ! Initialise the descriptor for each subproblem
      do i=1,nsubproblem
        call parlst_getvalue_string(rparlist, ssectionName, 'subproblem',&
            sproblemName, isubstring=i)
        call problem_initDescriptor(rproblemDescriptor%Rsubproblem(i),&
            rparlist, sproblemName)
      end do
      
      ! If this is a meta-problem we are done
      if (iand(rproblemDescriptor%cproblemtype,PROBTYPE_META) .ne. 0) return
    end if
    
    ! How many problem levels do we need?
    call parlst_getvalue_int(rparlist, ssectionName, 'nlmin',&
        rproblemDescriptor%nlmin)
    call parlst_getvalue_int(rparlist, ssectionName, 'nlmax',&
        rproblemDescriptor%nlmax)

    ! Next, we collect a data from the parameter list
    rproblemDescriptor%ndiscretisation = parlst_querysubstrings(rparlist,&
        ssectionName, 'discretisation')
    rproblemDescriptor%ncubatureInfo   = parlst_querysubstrings(rparlist,&
        ssectionName, 'cubatureinfo')
    rproblemDescriptor%nafcstab        = parlst_querysubstrings(rparlist,&
        ssectionName, 'afcstab')
    rproblemDescriptor%ngroupfemBlock  = parlst_querysubstrings(rparlist,&
        ssectionName, 'groupfemblock')
    rproblemDescriptor%nmatrixScalar   = parlst_querysubstrings(rparlist,&
        ssectionName, 'matrixscalar')
    rproblemDescriptor%nmatrixBlock    = parlst_querysubstrings(rparlist,&
        ssectionName, 'matrixblock')
    rproblemDescriptor%nvectorScalar   = parlst_querysubstrings(rparlist,&
        ssectionName, 'vectorscalar')
    rproblemDescriptor%nvectorBlock    = parlst_querysubstrings(rparlist,&
        ssectionName, 'vectorblock')
    
    ! Next, get triangulation data
    call parlst_getvalue_int(rparlist, ssectionName, 'ndimension',&
        rproblemDescriptor%ndimension)
    call parlst_getvalue_string(rparlist, ssectionName, 'trifile',&
        rproblemDescriptor%trifile, bdequote=.true.)
    call parlst_getvalue_string(rparlist, ssectionName, 'prmfile',&
        rproblemDescriptor%prmfile, '', bdequote=.true.)
    call parlst_getvalue_int(rparlist, ssectionName, 'iconvstrategy',&
        rproblemDescriptor%iconvStrategy, TRI_CONVERT_NONE)
    call parlst_getvalue_double(rparlist, ssectionName, 'ddisturbmesh',&
        rproblemDescriptor%ddisturbmesh, 0.0_DP)

#ifdef ENABLE_COPROCESSOR_SUPPORT
    ! Get number of CUDA streams
    call parlst_getvalue_int(rparlist, ssectionName, 'nstream',&
        rproblemDescriptor%nstream)
#endif

  end subroutine problem_initDescriptor

!*****************************************************************************

!<function>

  function problem_combineDescriptors(rproblemDescriptor1,&
                                      rproblemDescriptor2) result(rproblemDescriptor)

!<description>
    ! This function combines the given problem descriptors into a new
    ! one which comprises both descriptors.
!</description>
  
!<input>
    ! Problem descriptors to be combined
    type(t_problemDescriptor), intent(in) :: rproblemDescriptor1,rproblemDescriptor2
!</input>

!<result>
    ! Combined problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptor
!</result>
!</function>

    ! Spatial dimension
    if (rproblemDescriptor1%ndimension .ne. rproblemDescriptor2%ndimension) then
      call output_line('Spatial dimensions must be identical!',&
          OU_CLASS_WARNING,OU_MODE_STD,'problem_combineDescriptors')
      call sys_halt()
    else
      rproblemDescriptor%ndimension = rproblemDescriptor1%ndimension
    end if

    ! Triangulation
    if ((rproblemDescriptor1%trifile .ne. rproblemDescriptor2%trifile) .or.&
        (rproblemDescriptor1%prmfile .ne. rproblemDescriptor2%prmfile)) then
      call output_line('Triangulations must be identical!',&
          OU_CLASS_WARNING,OU_MODE_STD,'problem_combineDescriptors')
      call sys_halt()
    else
      rproblemDescriptor%trifile = rproblemDescriptor1%trifile
      rproblemDescriptor%prmfile = rproblemDescriptor1%prmfile
    end if
    
    ! Minimum/maximum problem level
    rproblemDescriptor%nlmax = max(rproblemDescriptor1%nlmax,&
                                   rproblemDescriptor2%nlmax)
    rproblemDescriptor%nlmin = min(rproblemDescriptor1%nlmin,&
                                   rproblemDescriptor2%nlmin)

    ! Discretisations, etc.
    rproblemDescriptor%ndiscretisation = max(rproblemDescriptor1%ndiscretisation,&
                                             rproblemDescriptor2%ndiscretisation)
    rproblemDescriptor%ncubatureInfo   = max(rproblemDescriptor1%ncubatureInfo,&
                                             rproblemDescriptor2%ncubatureInfo)
    rproblemDescriptor%nafcstab        = max(rproblemDescriptor1%nafcstab,&
                                             rproblemDescriptor2%nafcstab)
    rproblemDescriptor%ngroupfemBlock  = max(rproblemDescriptor1%ngroupfemBlock,&
                                             rproblemDescriptor2%ngroupfemBlock)
    rproblemDescriptor%nmatrixScalar   = max(rproblemDescriptor1%nmatrixScalar,&
                                             rproblemDescriptor2%nmatrixScalar)
    rproblemDescriptor%nmatrixBlock    = max(rproblemDescriptor1%nmatrixBlock,&
                                             rproblemDescriptor2%nmatrixBlock)
    rproblemDescriptor%nvectorScalar   = max(rproblemDescriptor1%nvectorScalar,&
                                             rproblemDescriptor2%nvectorScalar)
    rproblemDescriptor%nvectorBlock    = max(rproblemDescriptor1%nvectorBlock,&
                                             rproblemDescriptor2%nvectorBlock)


  end function problem_combineDescriptors

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_releaseDescriptor(rproblemDescriptor)

!<description>
    ! This subroutine releases the content of the problem descriptor
!</description>

!<inputoutput>
    ! Problem descriptor to be released
    type(t_problemDescriptor), intent(inout) :: rproblemDescriptor
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: i

    if (associated(rproblemDescriptor%Rsubproblem)) then
      do i = ubound(rproblemDescriptor%Rsubproblem,1),&
             lbound(rproblemDescriptor%Rsubproblem,1)
        call problem_releaseDescriptor(rproblemDescriptor%Rsubproblem(i))
      end do
      deallocate(rproblemDescriptor%Rsubproblem)
    end if

    call releaseContent(rproblemDescriptor)

  contains

    ! Here, the real info routines follow.

    !**************************************************************

    subroutine releaseContent(rproblemDescriptor)
      type(t_problemDescriptor), intent(out) :: rproblemDescriptor
      ! All data are released by INTENT(OUT)
    end subroutine releaseContent

  end subroutine problem_releaseDescriptor

  !*****************************************************************************

!<subroutine>

  subroutine problem_initDiscrAll(rproblemLevel, rparlist,&
      ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises all discretisation structures of the
    ! given problem level with the values provided by the parameter list.
    ! If the optional parameter p_rtasklist is given, then the task
    ! list which is created during the initialisation is returned.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! OPTIONAL: task list
    type(t_discrTask), intent(inout), pointer, optional :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_discrTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update all discretisations
    call problem_updateDiscrAll(rproblemLevel, rparlist,&
        ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearDiscrTL(p_rtasklistLocal)
    end if

  end subroutine problem_initDiscrAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_updateDiscrAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine generates the task list for all discretisation
    ! structures of the given problem level and updates all structures
    ! with the values provided by the parameter list.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! pointer to the task list
    type(t_discrTask), pointer :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: ssubsectionName
    integer :: i,ndiscretisation

    ! Consistency check
    ndiscretisation = parlst_querysubstrings(rparlist, ssectionName,&
                                             'discretisation')
    if (ndiscretisation .ne. size(rproblemLevel%Rdiscretisation)) then
      call output_line('Invalid number of discretisations',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscrAll')
      call sys_halt()
    end if

    ! Update all discretisations
    do i=1,ndiscretisation
      call parlst_getvalue_string(rparlist, ssectionName,&
          'discretisation', ssubsectionName, isubstring=i)
      call problem_updateDiscr(rproblemLevel,&
          rproblemLevel%Rdiscretisation(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do

  end subroutine problem_updateDiscrAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_initDiscr(rproblemLevel, rdiscretisation,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given discretisation structure
    ! on the given problem level with the values provided by the
    ! parameter list. If the optional parameter p_rtasklist is given,
    ! then the task list which is created during the initialisation is
    ! returned. Otherwise it is released on return.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which block
    ! discretisations need to be created and copied
    type(t_discrTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! block discretisation structure
    type(t_blockDiscretisation), intent(out), target :: rdiscretisation
!</output>
!</subroutine>

    ! local variables
    type(t_discrTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given discretisation
    call problem_updateDiscr(rproblemLevel, rdiscretisation,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearDiscrTL(p_rtasklistLocal)
    end if

  end subroutine problem_initDiscr

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateDiscr(rproblemLevel,&
      rdiscretisation, rparlist, ssectionName, p_rtasklist, rproblemTopLevel,&
      iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given discretisation structures on the given problem level with the
    ! values provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which block
    ! discretisations need to be created and copied
    type(t_discrTask), pointer :: p_rtasklist

    ! block discretisation structure
    type(t_blockDiscretisation), intent(inout), target :: rdiscretisation
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_discrTask), pointer :: rtask
    type(t_problem), pointer :: p_rproblemTopLevel
    type(t_problemLevel), pointer :: p_rproblemLevel
    character(len=PARLST_MLDATA) :: skeyword,sparameter,sproblem,stoken,svalue
    character(len=SYS_STRLEN) :: ssubsectionName
    character(len=SYS_STRLEN) :: saction,sperform,sperf
    integer :: idiscr,iistart,ilev,istart,itoken,ntoken,nblocks
    integer(i32) :: iperform

    ! Set pointer to problem structure
    p_rproblemTopLevel => rproblemLevel%p_rproblem
    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec

    ! Get type of action
    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
    call sys_toupper(saction)
    call sys_toupper(sperform)

    ! What action should be performed?
    if (trim(saction) .eq. 'NONE') then
      !-------------------------------------------------------------------------
      ! Do nothing
      
    elseif (trim(saction) .eq. 'CREATE') then

      !-------------------------------------------------------------------------
      ! Create block discretisation structure
      !
      ! SYNTAX: discretisation(n) =
      !           spatialdiscretisation = SpDiscr1
      !                                ...
      !           blockdiscretisation   = BlDiscr1
      !
      ! A block discretisation structure as the name suggests consists
      ! of multiple blocks of discretisation structures. Each block
      ! can be a spatial discretisation structure or a block
      ! discretisation structure. This hierarchic composition is
      ! traced down to spatial discretisation structures only.
      
      ! Create new task
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%ctask               =  PROBACTION_CREATE
      rtask%ssectionName        =  ssectionName
      rtask%p_rblockDiscrDest   => rdiscretisation
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelDest => rproblemLevel
      
      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)
          
          ! That`s it we do not have to create this discretisation
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
          call sys_halt()
        end if
      end do

      ! Create block discretisation structure
      if (iand(rtask%iperform, iperform) .ne. 0) then
        nblocks=0
        call countBlocks(rparlist, ssectionName, nblocks)
        call spdiscr_initBlockDiscr(rdiscretisation, nblocks,&
            rproblemLevel%rtriangulation, rproblemLevel%p_rproblem%rboundary)
        
        nblocks=0
        call initBlockDiscr(rparlist, ssectionName, rdiscretisation, nblocks,&
            rproblemLevel%rtriangulation, rproblemLevel%p_rproblem%rboundary)
      end if
      
    elseif (trim(saction) .eq. 'DUPLICATE') then

      !-------------------------------------------------------------------------
      ! Duplicate block discretisation structure
      !
      ! SYNTAX: discretisation = name,idiscretisation:#,ilev:#,...
      !                               ifirstblock:#,ilastblock:#
      !
      ! If ilev is not given then the level of the current problem
      ! level structure is adopted. If ifirstblock and/or ilastblock
      ! is not given then first and last block of the source block
      ! discretisation is used
 
      ! Find problem name from which to duplicate block discretisation
      call parlst_getvalue_string(rparlist, ssectionName, 'discretisation', sparameter)
      call sys_toupper(sparameter)

      ! Create new task for this discretisation structure
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%ctask = PROBACTION_DUPLICATE
      rtask%ssectionName = ssectionName

      ! Initialise optional parameters
      sproblem = trim(rproblemLevel%p_rproblem%cproblem)
      ilev     = rproblemLevel%ilev
      idiscr   = 1

      ! Get optional parameters if available
      call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
      
      do itoken = 1,ntoken
        call sys_getNextToken (sparameter, stoken, istart, ',', .false.)

        iistart = 1
        call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
        call sys_getNextToken (stoken, svalue, iistart, ":", .false.)

        if (trim(skeyword) .eq. 'PROBLEM') then
          sproblem = trim(svalue)
        elseif (trim(skeyword) .eq. 'ILEV') then
          read(svalue,'(I10)') ilev
        elseif (trim(skeyword) .eq. 'IFIRSTBLOCK') then
          read(svalue,'(I10)') rtask%ifirstBlock
        elseif (trim(skeyword) .eq. 'ILASTBLOCK') then
          read(svalue,'(I10)') rtask%ilastBlock
        elseif (trim(skeyword) .eq. 'IDISCR' .or.&
                trim(skeyword) .eq. 'IDISCRETISATION') then
          read(svalue,'(I10)') idiscr
        elseif (trim(skeyword) .eq. 'SHARE') then
          rtask%bshareStructure = (trim(svalue) .eq. 'YES')
        elseif (trim(skeyword) .eq. 'COPY') then
          rtask%bshareStructure = .not.(trim(svalue) .eq. 'YES')
        else
          call output_line('Unsupported keyword: '//trim(skeyword),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
          call sys_halt()
        end if
      end do

      ! Find problem level in global problem structure
      p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
      if (.not.associated(p_rproblemLevel)) then
        call output_line('Unable to find problem level in problem '//trim(sproblem),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
        call sys_halt()
      end if
      
      ! Get section name of discretisation
      call parlst_getvalue_string(rparlist, trim(sproblem),&
          'discretisation', ssubsectionName, isubstring=idiscr)
      
      ! Proceed with possibly prerequisite tasks
      call problem_updateDiscr(p_rproblemLevel,&
          p_rproblemLevel%Rdiscretisation(idiscr), rparlist,&
          ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
      
      ! Specify task for this discretisation structure
      rtask%p_rblockDiscrSrc    => p_rproblemLevel%Rdiscretisation(idiscr)
      rtask%p_rblockDiscrDest   => rdiscretisation
      rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelSrc  => p_rproblemLevel
      rtask%p_rproblemLevelDest => rproblemLevel

      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this discretisation
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
          call sys_halt()
        end if
      end do

      ! Duplicate block discretisation structure
      if ((iand(rtask%iperform, iperform) .ne. 0) .and.&
          (rtask%ifirstBlock .eq. 0) .and.&
          (rtask%ilastBlock .eq. 0)) then
        call spdiscr_duplicateBlockDiscr(rtask%p_rblockDiscrSrc,&
            rtask%p_rblockDiscrDest, rtask%bshareStructure)
      elseif (rtask%bshareStructure) then
        call spdiscr_deriveBlockDiscr(rtask%p_rblockDiscrSrc,&
            rtask%p_rblockDiscrDest, rtask%ifirstBlock, rtask%ilastBlock)
      else
        call output_line('Copying of a derived block discretisation is not supported',&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
        call sys_halt()
      end if
      
    else
      call output_line('Unsupported action: '//trim(saction),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
      call sys_halt()
    end if
    
  contains

    ! Here, some auxiliary routines follow

    !***************************************************************************
    ! Search for given task in the list of tasks. If the task is not
    ! present in the list, then it is appended.
    function appendToTaskList(rtasklist, rtask) result(bexists)

      type(t_discrTask), intent(inout), target :: rtasklist
      type(t_discrTask), intent(in), target :: rtask
      
      ! local variable
      type(t_discrTask), pointer :: p_rtask
      logical :: bexists      

      p_rtask => rtasklist
      do while(associated(p_rtask))
        if (associated(p_rtask%p_rblockDiscrDest,&
            rtask%p_rblockDiscrDest)) then
          bexists = .true.
          return
        end if
        if (.not.associated(p_rtask%p_rnextTask)) then
          p_rtask%p_rnextTask => rtask
          bexists = .false.
          return
        else
          p_rtask => p_rtask%p_rnextTask
        end if
      end do
      
    end function appendToTaskList
    
    !***************************************************************************
    ! Count total number of blocks in spatial discretisation
    recursive subroutine countBlocks(rparlist, ssectionName, nblocks)
      
      type(t_parlist), intent(in) :: rparlist
      character(len=*), intent(in) :: ssectionName
      integer, intent(inout) :: nblocks
      
      ! local variables
      character(len=PARLST_MLDATA) :: sparameter,stoken
      integer :: i,istart
      
      ! Loop over all substrings in "discretisation"
      do i = 0, parlst_querysubstrings(rparlist, ssectionName, "discretisation")
        call parlst_getvalue_string(rparlist, ssectionName, "discretisation",&
            sparameter, isubstring=i)
        
        istart = 1
        call sys_getNextToken (sparameter, stoken, istart, ":", .false.)
        call sys_toupper(stoken)
        
        ! Do we have a spatial- or blockdiscretisation structure?
        if (trim(adjustl(stoken)) .eq. "SPATIALDISCRETISATION") then
          nblocks = nblocks+1
          
        elseif (trim(adjustl(stoken)) .eq. "BLOCKDISCRETISATION") then
          call countBlocks(rparlist, trim(adjustl(sparameter(istart:))), nblocks)
          
        elseif(trim(adjustl(stoken)) .ne. "") then 
          call output_line('Unsupported parameter1: '//trim(adjustl(stoken)),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
          call sys_halt()
        end if
      end do
      
    end subroutine countBlocks

    !**************************************************************
    ! Initialise global block discretisation structure
    recursive subroutine initBlockDiscr(rparlist, ssectionName,&
        rdiscretisation, iblock, rtriangulation, rboundary)
      
      type(t_parlist), intent(in) :: rparlist
      character(len=*), intent(in) :: ssectionName
      type(t_triangulation), intent(in), target :: rtriangulation
      type(t_boundary), intent(in), target :: rboundary
      type(t_blockDiscretisation), intent(inout) :: rdiscretisation
      integer, dimension(:), allocatable :: Celement
      integer, intent(inout) :: iblock
      
      ! local variables
      character(len=PARLST_MLDATA) :: sparameter,stoken
      character(len=SYS_STRLEN) :: sspatialDiscrName,selemName
      integer :: i,ielement,istart,nelements
      
      ! Loop over all substrings in "discretisation"
      do i = 0, parlst_querysubstrings(rparlist, ssectionName, "discretisation")
        call parlst_getvalue_string(rparlist, ssectionName, "discretisation",&
            sparameter, isubstring=i)
        
        istart = 1
        call sys_getNextToken (sparameter, stoken, istart, ":", .false.)
        call sys_toupper(stoken)
        
        ! Do we have a spatial- or blockdiscretisation structure?
        if (trim(adjustl(stoken)) .eq. "SPATIALDISCRETISATION") then
          
          ! Increase block number
          iblock = iblock+1
          sspatialDiscrName = trim(adjustl(sparameter(istart:)))
          
          ! Allocate temporal memory
          nelements = max(1, parlst_querysubstrings(rparlist, sspatialDiscrName, 'celement'))
          allocate(Celement(nelements))
          
          ! Get IDs of element types
          do ielement = 1, nelements
            call parlst_getvalue_string(rparlist,&
                sspatialDiscrName, 'celement', selemName, isubstring=ielement)
            Celement(ielement) = elem_igetID(selemName)
          end do
          
          ! Get spatial dimension
          select case(rdiscretisation%ndimension)
          case (NDIM1D)
            call spdiscr_initDiscr_simple(rdiscretisation%RspatialDiscr(iblock),&
                Celement(1), rtriangulation, rboundary)
            
          case (NDIM2D)
            if (size(Celement) .eq. 1) then
              call spdiscr_initDiscr_simple(rdiscretisation%RspatialDiscr(iblock),&
                  Celement(1), rtriangulation, rboundary)
            else
              call spdiscr_initDiscr_triquad(&
                  rdiscretisation%RspatialDiscr(iblock), Celement(1), Celement(2),&
                  rtriangulation, rboundary)
            end if
            
          case (NDIM3D)
            call spdiscr_initDiscr_simple(rdiscretisation%RspatialDiscr(iblock),&
                Celement(1), rtriangulation, rboundary)
            
          case default
            call output_line('Invalid number of spatial dimensions',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
            call sys_halt()
          end select
          
          ! Deallocate temporal memory
          deallocate(Celement)

        elseif (trim(adjustl(stoken)) .eq. "BLOCKDISCRETISATION") then
          ! Proceed recursively for blockdiscretisation
          call initBlockDiscr(rparlist, trim(adjustl(sparameter(istart:))),&
              rdiscretisation, iblock, rtriangulation, rboundary)

        elseif(trim(adjustl(stoken)) .ne. "") then 
          call output_line('Unsupported parameter: '//trim(adjustl(stoken)),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateDiscr')
          call sys_halt()
        end if
      end do

    end subroutine initBlockDiscr
    
  end subroutine problem_updateDiscr
 
  !*****************************************************************************
  
!<subroutine>

  subroutine problem_clearDiscrTL(p_rtasklist)

!<description>
    ! This subroutine clears the discretisation task list.
!</description>

!<input>
    ! task list
    type(t_discrTask), pointer :: p_rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_discrTask), pointer :: p_rtask,p_rtaskSave

    p_rtask => p_rtasklist
    do while(associated(p_rtask))
      p_rtaskSave => p_rtask
      p_rtask     => p_rtask%p_rnextTask
      deallocate(p_rtaskSave)
    end do
    
  end subroutine problem_clearDiscrTL

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_infoDiscrTL(rtasklist)

!<description>
    ! This subroutine outputs information about the discretisation
    ! task list.
!</description>

!<input>
    ! task list
    type(t_discrTask), intent(in), target :: rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_discrTask), pointer :: p_rtask

    p_rtask => rtasklist
    do while(associated(p_rtask))

      select case(p_rtask%ctask)
      case (PROBACTION_NONE)
        call output_lbrk()
        call output_line ('DiscrTask: DO NOTHING')

      case (PROBACTION_CREATE)
        call output_lbrk()
        call output_line ('DiscrTask: CREATE')

      case (PROBACTION_DUPLICATE)
        call output_lbrk()
        call output_line ('DiscrTask: DUPLICATE')

      case default
        call output_line ('DiscrTask: UNSUPPORTED')
      end select

      call output_line ('----------------------')
      call output_line ('section name              : '//trim(p_rtask%ssectionName))
      call output_line ('destination problem       : '//trim(p_rtask%p_rproblemDest%cproblem))
      call output_line ('destination problem level : '//trim(sys_siL(p_rtask%p_rproblemLevelDest%ilev,15)))
      call output_line ('destination discretisation: '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rblockDiscrDest)))
      if (associated(p_rtask%p_rproblemSrc))&
      call output_line ('source problem            : '//trim(p_rtask%p_rproblemSrc%cproblem))
      if (associated(p_rtask%p_rproblemLevelSrc))&
      call output_line ('source problem level      : '//trim(sys_siL(p_rtask%p_rproblemLevelSrc%ilev,15)))
      call output_line ('source discretisation     : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rblockDiscrSrc)))
      call output_line ('iperform                  : '//trim(sys_siL(int(p_rtask%iperform),15)))
      call output_line ('shareStructure            : '//merge('TRUE ','FALSE',p_rtask%bshareStructure))
      call output_line ('ifirstBlock               : '//trim(sys_siL(p_rtask%ifirstBlock,15)))
      call output_line ('ilastBlock                : '//trim(sys_siL(p_rtask%ilastBlock,15)))

      p_rtask => p_rtask%p_rnextTask
    end do
    
  end subroutine problem_infoDiscrTL
  
  !*****************************************************************************

!<subroutine>

  subroutine problem_initCubInfoAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises all cubature info structures of the
    ! given problem level with the values provided by the parameter list.
    ! If the optional parameter p_rtasklist is given, then the task
    ! list which is created during the initialisation is returned.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! OPTIONAL: task list
    type(t_cubinfoTask), intent(inout), pointer, optional :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_cubinfoTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update all cubature info
    call problem_updateCubInfoAll(rproblemLevel,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearCubInfoTL(p_rtasklistLocal)
    end if

  end subroutine problem_initCubInfoAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_updateCubInfoAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine the task list for all cubature info structures of the
    ! given problem level with the values provided by the parameter list.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! pointer to the task list
    type(t_cubinfoTask), pointer :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: ssubsectionName
    integer :: i,ncubinfo

    ! Consistency check
    ncubinfo = parlst_querysubstrings(rparlist, ssectionName,&
                                      'cubatureinfo')
    if (ncubinfo .ne. size(rproblemLevel%RcubatureInfo)) then
      call output_line('Invalid number of cubature info structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfoAll')
      call sys_halt()
    end if

    ! Update all cubature info structures
    do i=1,ncubinfo
      call parlst_getvalue_string(rparlist, ssectionName,&
          'cubatureinfo', ssubsectionName, isubstring=i)
      call problem_updateCubInfo1(rproblemLevel,&
          rproblemLevel%Rcubatureinfo(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do

  end subroutine problem_updateCubInfoAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_initCubInfo(rproblemLevel, rcubatureinfo,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given cubature info structure on
    ! the given problem level with the values provided by the parameter
    ! list. If the optional parameter p_rtasklist is given, then the
    ! task list which is created during the initialisation is returned.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which cubature info
    ! structures need to be created and copied
    type(t_cubinfoTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! cubature info structure
    type(t_scalarCubatureInfo), intent(out), target :: rcubatureinfo
!</output>
!</subroutine>

    ! local variables
    type(t_cubinfoTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given cubature info structure
    call problem_updateCubInfo1(rproblemLevel, rcubatureinfo,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearCubInfoTL(p_rtasklistLocal)
    end if

  end subroutine problem_initCubInfo

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateCubInfo1(rproblemLevel,&
      rcubatureinfo, rparlist, ssectionName, p_rtasklist, rproblemTopLevel,&
      iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given cubature info structures on the given problem level with the
    ! values provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which cubature
    ! info structures need to be created and copied
    type(t_cubinfoTask), pointer :: p_rtasklist

    ! cubature info structure
    type(t_scalarCubatureInfo), intent(inout), target :: rcubatureinfo
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_cubinfoTask), pointer :: rtask
    type(t_problem), pointer :: p_rproblemTopLevel
    type(t_problemLevel), pointer :: p_rproblemLevel
    character(len=PARLST_MLDATA) :: skeyword,sparameter,scubtype
    character(len=PARLST_MLDATA) :: sproblem,stoken,svalue
    character(len=SYS_STRLEN) :: ssubsectionName
    character(len=SYS_STRLEN) :: saction,sperform,sperf
    integer :: icubinfo,iistart,ilev,istart,itoken,ntoken
    integer(i32) :: iperform

    ! Set pointer to problem structure
    p_rproblemTopLevel => rproblemLevel%p_rproblem
    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec

    ! Get type of action
    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
    call sys_toupper(saction)
    call sys_toupper(sperform)

    ! What action should be performed?
    if (trim(saction) .eq. 'NONE') then
      !-------------------------------------------------------------------------
      ! Do nothing

    elseif (trim(saction) .eq. 'CREATE') then

      !-------------------------------------------------------------------------
      ! Create cubature info structure
      !
      ! SYNTAX: cubatureinfo(n) =
      !           CubInfo1
      !             ...
      !           CubInfo2

      ! Create new task
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%ctask               =  PROBACTION_CREATE
      rtask%ssectionName        =  ssectionName
      rtask%p_rcubatureInfoDest => rcubatureinfo
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelDest => rproblemLevel
      
      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this cubature info structure
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo1')
          call sys_halt()
        end if
      end do

      ! Create cubature info structure
      if (iand(rtask%iperform, iperform) .ne. 0) then
        ! Get optional parameters
        call parlst_getvalue_int(rparlist, ssectionName,&
            'ccubType',  rtask%ccubType, CUB_UNDEFINED)
        if (rtask%ccubType .eq. CUB_UNDEFINED) then
          call parlst_getvalue_string(rparlist, ssectionName,&
              'scubType', scubType, 'AUTO')
          call sys_toupper(sparameter)
          rtask%ccubType = cub_igetID(scubtype)
        end if
        call parlst_getvalue_int(rparlist, ssectionName,&
            'nlevels', rtask%nlevels, 0)
        
        ! Get spatial discretisation
        call parlst_getvalue_string(rparlist, ssectionName,&
            'discretisation', sparameter)
        call sys_toupper(sparameter)
        rtask%p_rspatialDiscr => problem_getSpatialDiscr(&
            rtask%p_rproblemDest, rtask%p_rproblemLevelDest, sparameter)
        
        ! Do we have a spatial discretisation?
        if (.not.associated(rtask%p_rspatialDiscr)) then
          call output_line('Unable to find spatial discretisation',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo1')
          call sys_halt()
        end if
        
        ! Create cubature info structure
        call spdiscr_createDefCubStructure(rtask%p_rspatialDiscr,&
            rtask%p_rcubatureInfoDest, rtask%ccubType, rtask%nlevels)
      end if
            
    elseif (trim(saction) .eq. 'DUPLICATE') then

      !-------------------------------------------------------------------------
      ! Duplicate cubature info structure
      !
      ! SYNTAX: cubatureinfo = name,icubatureinfo:#,ilev:#,...
      !
      ! If ilev is not given then the level of the current problem
      ! level structure is adopted.

      ! Find problem name from which to duplicate cubature info
      call parlst_getvalue_string(rparlist, ssectionName, 'cubatureinfo', sparameter)
      call sys_toupper(sparameter)

      ! Create new task for this cubature info structure
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%ctask = PROBACTION_DUPLICATE
      rtask%ssectionName = ssectionName
      
      ! Initialise optional parameters
      sproblem = trim(rproblemLevel%p_rproblem%cproblem)
      ilev     = rproblemLevel%ilev
      icubinfo = 1
      
      ! Get optional parameters if available
      call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
      
      do itoken = 1,ntoken
        call sys_getNextToken (sparameter, stoken, istart, ',', .false.)

        iistart = 1
        call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
        call sys_getNextToken (stoken, svalue, iistart, ":", .false.)

        if (trim(skeyword) .eq. 'PROBLEM') then
          sproblem = trim(svalue)
        elseif (trim(skeyword) .eq. 'ILEV') then
          read(svalue,'(I10)') ilev
        elseif (trim(skeyword) .eq. 'ICUBINFO' .or.&
                trim(skeyword) .eq. 'ICUBATUREINRO') then
          read(svalue,'(I10)') icubinfo
        else
          call output_line('Unsupported keyword: '//trim(skeyword),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo1')
          call sys_halt()
        end if
      end do

      ! Find problem level in global problem structure
      p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
      if (.not.associated(p_rproblemLevel)) then
        call output_line('Unable to find problem level in problem '//trim(sproblem),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo1')
        call sys_halt()
      end if
          
      ! Get section name of cubature info structure
      call parlst_getvalue_string(rparlist, trim(sproblem),&
          'cubatureinfo', ssubsectionName, isubstring=icubinfo)
      
      ! Proceed with possibly prerequisite tasks
      call problem_updateCubInfo1(p_rproblemLevel,&
          p_rproblemLevel%Rcubatureinfo(icubinfo), rparlist,&
          ssubsectionName, p_rtasklist, p_rproblemTopLevel)
      
      ! Specify task for this cubature info structure
      rtask%p_rcubatureInfoSrc  => p_rproblemLevel%RcubatureInfo(icubinfo)
      rtask%p_rcubatureInfoDest => rcubatureinfo
      rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelSrc  => p_rproblemLevel
      rtask%p_rproblemLevelDest => rproblemLevel
      
      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this cubature info structure
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo1')
          call sys_halt()
        end if
      end do

      ! Duplicate cubature info structure
      if (iand(rtask%iperform, iperform) .ne. 0) then
        call output_line('Duplication of a cubature info structure is not supported ye',&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo1')
        call sys_halt()
      end if

    else
      call output_line('Unsupported action: '//trim(saction),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo1')
      call sys_halt()
    end if
    
  contains

    ! Here, some auxiliary routines follow

    !***************************************************************************
    ! Search for given task in the list of tasks. If the task is not
    ! present in the list, then it is appended.
    function appendToTaskList(rtasklist, rtask) result(bexists)

      type(t_cubinfoTask), intent(inout), target :: rtasklist
      type(t_cubinfoTask), intent(in), target :: rtask
      
      ! local variable
      type(t_cubinfoTask), pointer :: p_rtask
      logical :: bexists      

      p_rtask => rtasklist
      do while(associated(p_rtask))
        if (associated(p_rtask%p_rcubatureInfoDest,&
            rtask%p_rcubatureInfoDest)) then
          bexists = .true.
          return
        end if
        if (.not.associated(p_rtask%p_rnextTask)) then
          p_rtask%p_rnextTask => rtask
          bexists = .false.
          return
        else
          p_rtask => p_rtask%p_rnextTask
        end if
      end do
      
    end function appendToTaskList
    
  end subroutine problem_updateCubInfo1

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_updateCubInfo2(rtasklist, iperformSpec)

!<description>
    ! This subroutine updates the cubature info stuctures with the
    ! internal values assigned to the items of the task list.
!</description>

!<input>
    ! task list
    type(t_cubinfoTask), intent(in), target :: rtasklist

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>
!</subroutine>

    ! local variables
    type(t_cubinfoTask), pointer :: p_rtask
    integer :: iperform

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec
    
    ! Iterate over all tasks
    p_rtask => rtasklist
    do while (associated(p_rtask))

      ! Do we have to perform this task?
      if (iand(p_rtask%iperform, iperform) .eq. 0) then
        p_rtask => p_rtask%p_rnextTask
        cycle
      end if
      
      select case(p_rtask%ctask)
      case(PROBACTION_CREATE)
        
        !-----------------------------------------------------------------------
        ! Create cubature info structure

        call spdiscr_createDefCubStructure(p_rtask%p_rspatialDiscr,&
            p_rtask%p_rcubatureInfoDest, p_rtask%ccubType, p_rtask%nlevels)
        
      case(PROBACTION_DUPLICATE)

        !-----------------------------------------------------------------------
        ! Duplicate cubature info structure

        call output_line('Duplication of a cubature info structure is not supported ye',&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo2')
        call sys_halt()
        
      case default 
        call output_line('Unsupported action: '//sys_siL(p_rtask%ctask,3),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateCubInfo2')
        call sys_halt()
      end select
      
      ! Proceed with next task
      p_rtask => p_rtask%p_rnextTask
    end do
 
  end subroutine problem_updateCubInfo2

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_clearCubInfoTL(p_rtasklist)

!<description>
    ! This subroutine clears the cubature info task list.
!</description>

!<input>
    ! task list
    type(t_cubinfoTask), pointer :: p_rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_cubinfoTask), pointer :: p_rtask,p_rtaskSave

    p_rtask => p_rtasklist
    do while(associated(p_rtask))
      p_rtaskSave => p_rtask
      p_rtask     => p_rtask%p_rnextTask
      deallocate(p_rtaskSave)
    end do
    
  end subroutine problem_clearCubInfoTL

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_infoCubInfoTL(rtasklist)

!<description>
    ! This subroutine outputs information about the cubature info task
    ! list.
!</description>

!<input>
    ! task list
    type(t_cubinfoTask), intent(in), target :: rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_cubinfoTask), pointer :: p_rtask

    p_rtask => rtasklist
    do while(associated(p_rtask))
      
      select case(p_rtask%ctask)
      case (PROBACTION_NONE)
        call output_lbrk()
        call output_line ('CubInfoTask: DO NOTHING')

      case (PROBACTION_CREATE)
        call output_lbrk()
        call output_line ('CubInfoTask: CREATE')

      case (PROBACTION_DUPLICATE)
        call output_lbrk()
        call output_line ('CubInfoTask: DUPLICATE')

      case default
        call output_line ('CubInfoTask: UNSUPPORTED')
      end select

      call output_line ('------------------------')
      call output_line ('section name              : '//trim(p_rtask%ssectionName))
      call output_line ('destination problem       : '//trim(p_rtask%p_rproblemDest%cproblem))
      call output_line ('destination problem level : '//trim(sys_siL(p_rtask%p_rproblemLevelDest%ilev,15)))
      call output_line ('destination cubature info : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rcubatureInfoDest)))
      if (associated(p_rtask%p_rproblemSrc))&
      call output_line ('source problem            : '//trim(p_rtask%p_rproblemSrc%cproblem))
      if (associated(p_rtask%p_rproblemLevelSrc))&
      call output_line ('source problem level      : '//trim(sys_siL(p_rtask%p_rproblemLevelSrc%ilev,15)))
      call output_line ('source cubature info      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rcubatureInfoSrc)))
      call output_line ('spatial discretisation    : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rspatialDiscr)))
      call output_line ('iperform                  : '//trim(sys_siL(int(p_rtask%iperform),15)))
      call output_line ('cubature type             : '//trim(sys_siL(p_rtask%ccubType,15)))
      call output_line ('nlevels                   : '//trim(sys_siL(p_rtask%nlevels,15)))
      
      p_rtask => p_rtask%p_rnextTask
    end do
    
  end subroutine problem_infoCubInfoTL

  !*****************************************************************************

!<subroutine>

  subroutine problem_initMatrixAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises all scalar/block matrices of the
    ! given problem level with the values provided by the parameter
    ! list.  If the optional parameter p_rtasklist is given, then the
    ! task list which is created during the initialisation is returned.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! OPTIONAL: task list
    type(t_matrixTask), intent(inout), pointer, optional :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_matrixTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update all scalar/block matrices
    call problem_updateMatrixAll(rproblemLevel,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearMatrixTL(p_rtasklistLocal)
    end if

  end subroutine problem_initMatrixAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_updateMatrixAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine the task list for all scalar/block matrices of the
    ! given problem level with the values provided by the parameter list.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! pointer to the task list
    type(t_matrixTask), pointer :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: ssubsectionName
    integer :: i,nmatrixScalar,nmatrixBlock

    ! Consistency check
    nmatrixScalar = parlst_querysubstrings(rparlist, ssectionName,&
                                           'matrixscalar')
    if ((nmatrixScalar .gt. 0) .and.&
        (nmatrixScalar .ne. size(rproblemLevel%RmatrixScalar))) then
      call output_line('Invalid number of scalar matrices',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixAll')
      call sys_halt()
    end if

    nmatrixBlock = parlst_querysubstrings(rparlist, ssectionName,&
                                          'matrixblock')
    if ((nmatrixBlock .gt. 0) .and.&
        (nmatrixBlock .ne. size(rproblemLevel%RmatrixBlock))) then
      call output_line('Invalid number of block matrices',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixAll')
      call sys_halt()
    end if
    
    ! Update all scalar matrices
    do i=1,nmatrixScalar
      call parlst_getvalue_string(rparlist, ssectionName,&
          'matrixscalar', ssubsectionName, isubstring=i)
      call problem_updateMatrixScalar(rproblemLevel,&
          rproblemLevel%RmatrixScalar(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do
    
    ! Update all block matrices
    do i=1,nmatrixBlock
      call parlst_getvalue_string(rparlist, ssectionName,&
          'matrixblock', ssubsectionName, isubstring=i)
      call problem_updateMatrixBlock(rproblemLevel,&
          rproblemLevel%RmatrixBlock(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do

  end subroutine problem_updateMatrixAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_initMatrixScalar(rproblemLevel, rmatrix,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given scalar matrix on the given
    ! problem level with the values provided by the parameter list. If
    ! the optional parameter p_rtasklist is given, then the task list
    ! which is created during the initialisation is returned.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar/block
    ! matrices need to be created and copied
    type(t_matrixTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! scalar matrix
    type(t_matrixScalar), intent(out), target :: rmatrix
!</output>
!</subroutine>

    ! local variables
    type(t_matrixTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given scalar matrix
    call problem_updateMatrixScalar(rproblemLevel, rmatrix,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearMatrixTL(p_rtasklistLocal)
    end if

  end subroutine problem_initMatrixScalar

  !*****************************************************************************

!<subroutine>

  subroutine problem_initMatrixBlock(rproblemLevel, rmatrix,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given block matrix on the given
    ! problem level with the values provided by the parameter list. If
    ! the optional parameter p_rtasklist is given, then the task list
    ! which is created during the initialisation is returned.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar/block
    ! matrices need to be created and copied
    type(t_matrixTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! block matrix
    type(t_matrixBlock), intent(out), target :: rmatrix
!</output>
!</subroutine>

    ! local variables
    type(t_matrixTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given block matrix
    call problem_updateMatrixBlock(rproblemLevel, rmatrix,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearMatrixTL(p_rtasklistLocal)
    end if

  end subroutine problem_initMatrixBlock

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateMatrixScalar(rproblemLevel,&
      rmatrix, rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given scalar matrix on the given problem level with the values
    ! provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar matrices
    ! need to be created and copied
    type(t_matrixTask), pointer :: p_rtasklist

    ! scalar matrix
    type(t_matrixScalar), intent(inout), target :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock), pointer :: p_rmatrixBlock
    type(t_matrixTask), pointer :: rtask
    type(t_problem), pointer :: p_rproblemTopLevel
    type(t_problemLevel), pointer :: p_rproblemLevel
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    character(len=PARLST_MLDATA) :: skeyword,sparameter,sproblem,stoken,svalue
    character(len=SYS_STRLEN) :: ssubsectionName
    character(len=SYS_STRLEN) :: saction,sperform,sperf
    integer :: iblockcol,iblockrow,iistart,ilev,imatrix,istart,itoken,ntoken,imethod,nmethod
    integer(i32) :: iperform

    ! Set pointer to problem structure
    p_rproblemTopLevel => rproblemLevel%p_rproblem
    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec

    ! Get type of action
    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
    call sys_toupper(saction)
    call sys_toupper(sperform)

    !---------------------------------------------------------------------------
    ! First part: generate structure
    !---------------------------------------------------------------------------

    ! What action should be performed?
    if (trim(saction) .eq. 'NONE') then
      !-------------------------------------------------------------------------
      ! Do nothing

    elseif (trim(saction) .eq. 'CREATE') then

      !-------------------------------------------------------------------------
      ! Create scalar matrix
      !
      ! SYNTAX: matrixscalar(n) =
      !           MatrixScalar1
      !                ...
      !           MatrixScalar2

      ! Create new task
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%bisMatrixScalar     = .true.
      rtask%ctask               =  PROBACTION_CREATE
      rtask%ssectionName        =  ssectionName
      rtask%p_rmatrixScalarDest => rmatrix
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelDest => rproblemLevel

      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this scalar matrix
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
          call sys_halt()
        end if
      end do
      
      ! Get optional parameters
      call  parlst_getvalue_int(rparlist, ssectionName,&
          'cmatrixformat', rtask%cmatrixFormat, LSYSSC_MATRIXUNDEFINED)
      if (rtask%cmatrixFormat .eq. LSYSSC_MATRIXUNDEFINED) then
        call parlst_getvalue_string(rparlist, ssectionName,&
            'smatrixformat', sparameter, 'MATRIX9')
        call sys_toupper(sparameter)
        rtask%cmatrixFormat = get_matrixFormat(trim(sparameter))
      end if
      
      call parlst_getvalue_int(rparlist, ssectionName,&
          'cinterleavematrixformat', rtask%cinterleavematrixFormat, LSYSSC_MATRIXUNDEFINED)
      if (rtask%cinterleavematrixFormat .eq. LSYSSC_MATRIXUNDEFINED) then
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sinterleavematrixformat', sparameter, 'MATRIXUNDEFINED')
        call sys_toupper(sparameter)
        rtask%cinterleavematrixFormat = get_interleavematrixFormat(trim(sparameter))
      end if

      call parlst_getvalue_int(rparlist, ssectionName,&
          'cdatatype', rtask%cdataType, ST_NOHANDLE)
      if (rtask%cdataType .eq. ST_NOHANDLE) then
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sdatatype', sparameter, 'DOUBLE')
        call sys_toupper(sparameter)
        rtask%cdataType = get_dataType(trim(sparameter))
      end if

      call parlst_getvalue_int(rparlist, ssectionName,&
          'nvar', rtask%nvar, 1)
      call parlst_getvalue_double(rparlist, ssectionName,&
          'dscalefactor', rtask%dscaleFactor, 1.0_DP)

      if ((rtask%cmatrixFormat .ne. LSYSSC_MATRIX7INTL) .and.&
          (rtask%cmatrixFormat .ne. LSYSSC_MATRIX9INTL)) then
        rtask%nvar = 1
      end if
      
      ! Get spatial discretisation for test and trial space
      call parlst_getvalue_string(rparlist, ssectionName,&
          'discretisation', sparameter, '')
      if (trim(sparameter) .eq. '') then
        ! Get spatial discretisation for test space
        call parlst_getvalue_string(rparlist, ssectionName,&
            'discretisationtest', sparameter)
        call sys_toupper(sparameter)
        rtask%p_rspatialDiscrTest => problem_getSpatialDiscr(&
            rproblemLevel%p_rproblem, rproblemLevel, sparameter)
        ! Get spatial discretisation for trial space
        call parlst_getvalue_string(rparlist, ssectionName,&
            'discretisationtrial', sparameter)
        call sys_toupper(sparameter)
        rtask%p_rspatialDiscrTrial => problem_getSpatialDiscr(&
            rproblemLevel%p_rproblem, rproblemLevel, sparameter)
      else
        call sys_toupper(sparameter)
        p_rspatialDiscr => problem_getSpatialDiscr(&
            rproblemLevel%p_rproblem, rproblemLevel, sparameter)
        rtask%p_rspatialDiscrTest => p_rspatialDiscr
        rtask%p_rspatialDiscrTrial => p_rspatialDiscr
      end if

      ! Should we perform this task now?
      if (iand(rtask%iperform, iperform) .ne. 0) then
        ! Do we have spatial discretisations?
        if (.not.associated(rtask%p_rspatialDiscrTest) .or.&
            .not.associated(rtask%p_rspatialDiscrTrial)) then
          call output_line('Unable to find spatial discretisation',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
          call sys_halt()
        end if
        
        ! Clear preexisting scalar matrix
        call lsyssc_releaseMatrix(rmatrix)
        
        ! Create matrix structure and set data type
        if (associated(rtask%p_rspatialDiscrTest,rtask%p_rspatialDiscrTrial)) then
          call bilf_createMatrixStructure(rtask%p_rspatialDiscrTrial,&
              rtask%cmatrixFormat, rmatrix)
        else
          call bilf_createMatrixStructure(rtask%p_rspatialDiscrTrial,&
              rtask%cmatrixFormat, rmatrix, rtask%p_rspatialDiscrTest)
        end if
        call lsyssc_setDataTypeMatrix (rmatrix, rtask%cdataType)
        rmatrix%dscaleFactor = rtask%dscaleFactor
        if ((rtask%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) .or.&
            (rtask%cmatrixFormat .eq. LSYSSC_MATRIX9INTL)) then
          rmatrix%nvar = rtask%nvar
        end if
      end if

    elseif (trim(saction) .eq. 'DUPLICATE') then

      !-------------------------------------------------------------------------
      ! Duplicate scalar matrix or submatrix of a block matrix
      !
      ! SYNTAX: matrixscalar = name,imatrix:#,ilev:#,...
      !     OR  matrixblock  = name,imatrix:#,ilev:#,iblockrow:#,iblockcol:#
      !
      ! If ilev is not given then the level of the current problem
      ! level structure is adopted. If iblockrow and/or iblockcol are
      ! not given then standard value 1 is used in both cases

      ! Find problem name from which to duplicate scalar matrix (if any)
      call parlst_getvalue_string(rparlist, ssectionName, 'matrixscalar', sparameter, '')
      call sys_toupper(sparameter)

      if (trim(sparameter) .ne. '') then
        
        !--- Scalar matrix case ------------------------------------------------

        ! Create new task for this scalar matrix
        nullify(rtask)
        allocate(rtask)
        nullify(rtask%p_rnextTask)
        
        ! Initialise optional parameters
        sproblem = trim(rproblemLevel%p_rproblem%cproblem)
        ilev     = rproblemLevel%ilev
        imatrix  = 1
        
        ! Get optional parameters if available
        call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
        
        do itoken = 1,ntoken
          call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
          
          iistart = 1
          call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
          call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
          
          if (trim(skeyword) .eq. 'PROBLEM') then
            sproblem = trim(svalue)
          elseif (trim(skeyword) .eq. 'ILEV') then
            read(svalue,'(I10)') ilev
          elseif (trim(skeyword) .eq. 'IMATRIX') then
            read(svalue,'(I10)') imatrix
          else
            call output_line('Unsupported keyword: '//trim(skeyword),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
            call sys_halt()
          end if
        end do
        
        ! Find problem level in global problem structure
        p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
        if (.not.associated(p_rproblemLevel)) then
          call output_line('Unable to find problem level in problem '//trim(sproblem),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
          call sys_halt()
        end if
        
        ! Get section name of scalar matrix
        call parlst_getvalue_string(rparlist, trim(sproblem),&
            'matrixscalar', ssubsectionName, isubstring=imatrix)

        ! Proceed with possibly prerequisite tasks
        call problem_updateMatrixScalar(p_rproblemLevel,&
            p_rproblemLevel%RmatrixScalar(imatrix), rparlist,&
            ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
        
        ! Specify new task for this scalar matrix
        rtask%bisMatrixScalar     = .true.
        rtask%ctask               =  PROBACTION_DUPLICATE
        rtask%ssectionName        =  ssectionName
        rtask%p_rmatrixScalarSrc  => p_rproblemLevel%RmatrixScalar(imatrix)
        rtask%p_rmatrixScalarDest => rmatrix
        rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
        rtask%p_rproblemDest      => rproblemLevel%p_rproblem
        rtask%p_rproblemLevelSrc  => p_rproblemLevel
        rtask%p_rproblemLevelDest => rproblemLevel

      else

        !--- Block matrix case -------------------------------------------------

        ! Find problem name from which to duplicate entry of block matrix (if any)
        call parlst_getvalue_string(rparlist, ssectionName, 'matrixblock', sparameter, '')
        call sys_toupper(sparameter)
      
        if (trim(sparameter) .ne. '') then
          
          ! Create new task for this scalar matrix
          nullify(rtask)
          allocate(rtask)
          nullify(rtask%p_rnextTask)
          
          ! Initialise optional parameters
          sproblem = trim(rproblemLevel%p_rproblem%cproblem)
          ilev     = rproblemLevel%ilev
          imatrix  = 1
          iblockrow = 1
          iblockcol = 1
          
          ! Get optional parameters if available
          call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
          
          do itoken = 1,ntoken
            call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
            
            iistart = 1
            call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
            call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
            
            if (trim(skeyword) .eq. 'PROBLEM') then
              sproblem = trim(svalue)
            elseif (trim(skeyword) .eq. 'ILEV') then
              read(svalue,'(I10)') ilev
            elseif (trim(skeyword) .eq. 'IMATRIX') then
              read(svalue,'(I10)') imatrix
            elseif (trim(skeyword) .eq. 'IBLOCKROW') then
              read(svalue,'(I10)') iblockrow
            elseif (trim(skeyword) .eq. 'IBLOCKCOL') then
              read(svalue,'(I10)') iblockcol
            else
              call output_line('Unsupported keyword: '//trim(skeyword),&
                  OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
              call sys_halt()
            end if
          end do
          
          ! Find problem level in global problem structure
          p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
          if (.not.associated(p_rproblemLevel)) then
            call output_line('Unable to find problem level in problem '//trim(sproblem),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
            call sys_halt()
          end if
          
          ! Get section name of block matrix
          call parlst_getvalue_string(rparlist, trim(sproblem),&
              'matrixblock', ssubsectionName, isubstring=imatrix)

          ! Proceed with possibly prerequisite tasks
          call problem_updateMatrixBlock(p_rproblemLevel,&
              p_rproblemLevel%RmatrixBlock(imatrix), rparlist,&
              ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)

          ! Check if scalar submatrix is available
          p_rmatrixBlock => p_rproblemLevel%RmatrixBlock(imatrix)
          if ((size(p_rmatrixBlock%RmatrixBlock,1) .lt. iblockrow) .or.&
              (size(p_rmatrixBlock%RmatrixBlock,2) .lt. iblockcol)) then
            call output_line('Scalar submatrix of block matrix is not available',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
            call sys_halt()
          end if

          ! Specify new task for this scalar matrix
          rtask%bisMatrixScalar     = .true.
          rtask%ctask               =  PROBACTION_DUPLICATE
          rtask%ssectionName        =  ssectionName
          rtask%p_rmatrixScalarSrc  => p_rmatrixBlock%RmatrixBlock(iblockrow,iblockcol)
          rtask%p_rmatrixScalarDest => rmatrix
          rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
          rtask%p_rproblemDest      => rproblemLevel%p_rproblem
          rtask%p_rproblemLevelSrc  => p_rproblemLevel
          rtask%p_rproblemLevelDest => rproblemLevel

        else
          call output_line('Either matrixscalar of matrixblock must be present',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
          call sys_halt()
        end if
      end if

      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this scalar matrix
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
          call sys_halt()
        end if
      end do
          
      ! Get optional parameters
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupstructure', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupStructure = get_dupType(sparameter)
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupcontent', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupContent = get_dupType(sparameter)

      ! Duplicate scalar matrix
      if (iand(rtask%iperform, iperform) .ne. 0) then
        call lsyssc_duplicateMatrix(rtask%p_rmatrixScalarSrc,&
            rmatrix, rtask%cdupStructure, rtask%cdupContent)
      end if
      
    else
      call output_line('Unsupported action: '//trim(saction),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
      call sys_halt()
    end if
    
    !---------------------------------------------------------------------------
    ! Second part: generate content
    !---------------------------------------------------------------------------

    nmethod = parlst_querysubstrings(rparlist, ssectionName, 'smethod')
    allocate(rtask%Cmethod(0:nmethod))

    ! What creation methods are adopted?
    do imethod = 0,nmethod
      
      call parlst_getvalue_string(rparlist, ssectionName,&
          'smethod', sparameter, '', isubstring=imethod)
      call sys_toupper(sparameter)
      
      if (trim(sparameter) .eq. 'VIRTUAL') then
        
        !-----------------------------------------------------------------------
        ! Create virtual matrix, aka, do nothing
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_VIRTUAL
        
      elseif (trim(sparameter) .eq. 'EMPTY') then
        
        !-----------------------------------------------------------------------
        ! Create empty matrix
        
        call lsyssc_allocEmptyMatrix(rmatrix, LSYSSC_SETM_ZERO)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_EMPTY
        
      elseif (trim(sparameter) .eq. 'IDENTITY') then
        
        !-----------------------------------------------------------------------
        ! Create identity matrix
        
        call lsyssc_initialiseIdentityMatrix(rmatrix)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_IDENTITY
        
      elseif (trim(sparameter) .eq. 'BILF') then
        
        !-----------------------------------------------------------------------
        ! Create matrix by bilinearform evaluation
        
        call assembleScalarMatrixBilf(rparlist, rtask)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_BILF
        
      elseif (trim(sparameter) .eq. 'LUMP_STD') then
        
        !-----------------------------------------------------------------------
        ! Create matrix by standard lumping
        
        call lsyssc_lumpMatrixScalar(rmatrix, LSYSSC_LUMP_STD)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_LUMP_STD
        
      elseif (trim(sparameter) .eq. 'LUMP_DIAG') then
        
        !-----------------------------------------------------------------------
        ! Create matrix by diagonal lumping
        
        call lsyssc_lumpMatrixScalar(rmatrix, LSYSSC_LUMP_DIAG)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_LUMP_DIAG
        
      elseif (trim(sparameter) .eq. 'EXTENDSPARSITY') then
        
        !-----------------------------------------------------------------------
        ! Create matrix with extended sparsity pattern
        
        call afcstab_genExtSparsity(rmatrix)

      else
        call output_line('Unsupported method: '//trim(sparameter),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixScalar')
        call sys_halt()
      end if
    end do

  contains

    ! Here, some auxiliary routines follow

    !***************************************************************************
    ! Search for given task in the list of tasks. If the task is not
    ! present in the list, then it is appended.
    function appendToTaskList(rtasklist, rtask) result(bexists)

      type(t_matrixTask), intent(inout), target :: rtasklist
      type(t_matrixTask), intent(in), target :: rtask
      
      ! local variable
      type(t_matrixTask), pointer :: p_rtask
      logical :: bexists      

      p_rtask => rtasklist
      do while(associated(p_rtask))
        if (associated(p_rtask%p_rmatrixScalarDest,&
            rtask%p_rmatrixScalarDest)) then
          bexists = .true.
          return
        end if
        if (.not.associated(p_rtask%p_rnextTask)) then
          p_rtask%p_rnextTask => rtask
          bexists = .false.
          return
        else
          p_rtask => p_rtask%p_rnextTask
        end if
      end do
      
    end function appendToTaskList

    !**************************************************************
    ! Return the numeric value of the data type
    function get_dataType(sdataType) result(cdataType)

      character(len=*), intent(in) :: sdataType
      integer :: cdataType

      if (sdataType .eq. 'QUAD') then
        cdataType = ST_QUAD
      elseif (sdataType .eq. 'DOUBLE') then
        cdataType = ST_DOUBLE
      elseif (sdataType .eq. 'SINGLE') then
        cdataType = ST_SINGLE
      else
        cdataType = ST_NOHANDLE
      end if

    end function get_dataType

    !**************************************************************
    ! Return the numeric value of the matrix format
    function get_matrixFormat(smatrixFormat) result(cmatrixFormat)

      character(len=*), intent(in) :: smatrixFormat
      integer :: cmatrixFormat

      if (smatrixFormat .eq. 'MATRIX1') then
        cmatrixFormat = LSYSSC_MATRIX1
      elseif (smatrixFormat .eq. 'MATRIX7') then
        cmatrixFormat = LSYSSC_MATRIX7
      elseif (smatrixFormat .eq. 'MATRIX7INTL') then
        cmatrixFormat = LSYSSC_MATRIX7INTL
      elseif (smatrixFormat .eq. 'MATRIX9') then
        cmatrixFormat = LSYSSC_MATRIX9
      elseif (smatrixFormat .eq. 'MATRIX9INTL') then
        cmatrixFormat = LSYSSC_MATRIX9INTL
      elseif (smatrixFormat .eq. 'MATRIXD') then
        cmatrixFormat = LSYSSC_MATRIXD
      else
        cmatrixFormat = LSYSSC_MATRIXUNDEFINED
      end if

    end function get_matrixFormat
    
    !**************************************************************
    ! Return the numeric value of the interleaved matrix format
    function get_interleavematrixFormat(sinterleavematrixFormat)&
        result(cinterleavematrixFormat)
      
      character(len=*), intent(in) :: sinterleavematrixFormat
      integer :: cinterleavematrixFormat
      
      if (sinterleavematrixFormat .eq. 'MATRIX1') then
        cinterleavematrixFormat = LSYSSC_MATRIX1
      elseif (sinterleavematrixFormat .eq. 'MATRIXD') then
        cinterleavematrixFormat = LSYSSC_MATRIXD
      else
        cinterleavematrixFormat = LSYSSC_MATRIXUNDEFINED
      end if
      
    end function get_interleavematrixFormat

    !**************************************************************
    ! Return the numeric value of the duplication flad
    function get_dupType(sdupType) result(cdupType)

      character(len=*), intent(in) :: sdupType
      integer :: cdupType

      if (sdupType .eq. 'IGNORE') then
        cdupType = LSYSSC_DUP_IGNORE
      elseif (sdupType .eq. 'REMOVE') then
        cdupType = LSYSSC_DUP_REMOVE
      elseif (sdupType .eq. 'DISMISS') then
        cdupType = LSYSSC_DUP_DISMISS
      elseif (sdupType .eq. 'SHARE') then
        cdupType = LSYSSC_DUP_SHARE
      elseif (sdupType .eq. 'COPY') then
        cdupType = LSYSSC_DUP_COPY
      elseif (sdupType .eq. 'COPYOVERWRITE') then
        cdupType = LSYSSC_DUP_COPYOVERWRITE
      elseif (sdupType .eq. 'ASIS') then
        cdupType = LSYSSC_DUP_ASIS
      elseif (sdupType .eq. 'EMPTY') then
        cdupType = LSYSSC_DUP_EMPTY
      elseif (sdupType .eq. 'TEMPLATE') then
        cdupType = LSYSSC_DUP_TEMPLATE
      else
        cdupType = 0
      end if

    end function get_dupType
    
    !**************************************************************
    ! Assemble the scalar matrix from bilinear form evaluation
    subroutine assembleScalarMatrixBilf(rparlist, rtask)

      type(t_parlist), intent(in) :: rparlist
      type(t_matrixTask), intent(inout) :: rtask

      ! local variables
      type(t_collection) :: rcollection
      type(t_scalarCubatureInfo), pointer :: p_rcubInfo
      character(len=PARLST_MLDATA) :: sparameter,stoken
      integer :: i,j,itermCount,istart

      ! Get cubature info structure
      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
          'cubatureinfo', sparameter)
      call sys_toupper(sparameter)

      p_rcubInfo => problem_getCubInfo(rtask%p_rproblemDest,&
          rtask%p_rproblemLevelDest, sparameter)
      
      ! Update internal data of the task item
      rtask%p_rcubatureInfo => p_rcubInfo
      allocate(rtask%p_rform)

      ! Setup trial functions
      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
          'strialfunction', sparameter)
      call sys_countTokens(sparameter, itermCount, ',', .false.)
      rtask%p_rform%itermCount = itermCount
      
      istart = 1
      do i = 1,itermCount
        call sys_getNextToken(sparameter, stoken, istart, ',', .false.)
        rtask%p_rform%Idescriptors(1,i) = der_igetID(stoken)
      end do
      
      ! Setup test functions
      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
          'stestfunction', sparameter)
      call sys_countTokens(sparameter, itermCount, ',', .false.)
      rtask%p_rform%itermCount = min(rtask%p_rform%itermCount,itermCount)
      
      istart = 1
      do i = 1,itermCount
        call sys_getNextToken(sparameter, stoken, istart, ',', .false.)
        rtask%p_rform%Idescriptors(2,i) = der_igetID(stoken)
      end do

      ! Setup coefficients
      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
          'scoefficient', sparameter, '')
      call sys_countTokens(sparameter, itermCount, ',', .false.)
      
      if (itermCount .gt. 0) then
        ! Use individual coefficients for each component
        rtask%p_rform%itermCount = min(rtask%p_rform%itermCount,itermCount)

        istart = 1
        do i=1,itermCount
          call sys_getNextToken(sparameter, stoken, istart, ',', .false.)
          if (sys_isNumeric(stoken)) then
            read(stoken,*) rtask%p_rform%Dcoefficients(i)
          else
            ! Not all coefficients are constant. Thus, we create a
            ! function parser object and perform bilinear form
            ! evaluation based on a generic callback function
            rtask%p_rform%ballCoeffConstant = .false.
            rtask%p_rform%BconstantCoeff(1:itermCount) = .false.
            
            allocate (rtask%p_rfparser)
            call fparser_create(rtask%p_rfparser, itermCount)
            
            ! Fill function parser with data
            istart = 1
            do j = 1,itermCount
              call sys_getNextToken(sparameter, stoken, istart, ',', .false.)
              call fparser_parseFunction(rtask%p_rfparser, j, trim(stoken), (/'x','y','z'/))
            end do
            
            ! That`s it
            exit
          end if
        end do
        
      else
        ! Use unity as constand coefficient
        rtask%p_rform%Dcoefficients(1:rtask%p_rform%itermCount) = 1.0_DP
      end if
          
      ! Assemble matrix data from bilinear form
      if (rtask%p_rform%ballCoeffConstant) then
        if (associated(p_rcubInfo)) then
          call bilf_buildMatrixScalar(rtask%p_rform, .true.,&
              rtask%p_rmatrixScalarDest, p_rcubInfo)
        else
          call bilf_buildMatrixScalar(rtask%p_rform, .true.,&
              rtask%p_rmatrixScalarDest)
        end if
      else
        rcollection%p_rfparserQuickAccess1 => rtask%p_rfparser
        if (associated(p_rcubInfo)) then
          call bilf_buildMatrixScalar(rtask%p_rform, .true.,&
              rtask%p_rmatrixScalarDest, p_rcubInfo,&
              problem_buildMatrixSc_fparser_sim, rcollection)
        else
          call bilf_buildMatrixScalar(rtask%p_rform, .true.,&
              rtask%p_rmatrixScalarDest,&
              problem_buildMatrixSc_fparser_sim, rcollection)
        end if
      end if
    
    end subroutine assembleScalarMatrixBilf

  end subroutine problem_updateMatrixScalar

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateMatrixBlock(rproblemLevel,&
      rmatrix, rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given block matrix on the given problem level with the values
    ! provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar matrices
    ! need to be created and copied
    type(t_matrixTask), pointer :: p_rtasklist

    ! block matrix
    type(t_matrixBlock), intent(inout), target :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock) :: rmatrixTemp 
    type(t_matrixTask), pointer :: rtask
    type(t_problem), pointer :: p_rproblemTopLevel
    type(t_problemLevel), pointer :: p_rproblemLevel
    type(t_blockDiscretisation), pointer :: p_rblockDiscr
    character(len=PARLST_MLDATA) :: skeyword,sparameter,sproblem,stoken,svalue
    character(len=SYS_STRLEN) :: ssubsectionName
    character(len=SYS_STRLEN) :: saction,sperform,sperf
    integer :: iistart,ilev,imatrix,istart,itoken,ntoken,i,j,imethod,nmethod
    integer(i32) :: iperform

    ! Set pointer to problem structure
    p_rproblemTopLevel => rproblemLevel%p_rproblem
    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec

    ! Get type of action
    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
    call sys_toupper(saction)
    call sys_toupper(sperform)

    ! What action should be performed?
    if (trim(saction) .eq. 'NONE') then
      !-------------------------------------------------------------------------
      ! Do nothing
      
    elseif (trim(saction) .eq. 'CREATE') then

      !-------------------------------------------------------------------------
      ! Create block matrix
      !
      ! SYNTAX: matrixblock(n) =
      !           MatrixBlock1
      !                ...
      !           MatrixBlock2

      ! Create new task
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%bisMatrixScalar     = .false.
      rtask%ctask               =  PROBACTION_CREATE
      rtask%ssectionName        =  ssectionName
      rtask%p_rmatrixBlockDest  => rmatrix
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelDest => rproblemLevel
      
      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this block matrix
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
          call sys_halt()
        end if
      end do

      ! Get block discretisation for test and trial space
      call parlst_getvalue_string(rparlist, ssectionName,&
          'discretisation', sparameter, '')
      if (trim(sparameter) .eq. '') then
        ! Get block discretisation for test space
        call parlst_getvalue_string(rparlist, ssectionName,&
            'discretisationtest', sparameter)
        call sys_toupper(sparameter)
        rtask%p_rblockDiscrTest => problem_getBlockDiscr(rproblemLevel%p_rproblem,&
            rproblemLevel, sparameter)
        ! Get block discretisation for trial space
        call parlst_getvalue_string(rparlist, ssectionName,&
            'discretisationtrial', sparameter)
        call sys_toupper(sparameter)
        rtask%p_rblockDiscrTrial => problem_getBlockDiscr(rproblemLevel%p_rproblem,&
            rproblemLevel, sparameter)
      else
        call sys_toupper(sparameter)
        p_rblockDiscr => problem_getBlockDiscr(rproblemLevel%p_rproblem,&
            rproblemLevel, sparameter)
        rtask%p_rblockDiscrTest => p_rblockDiscr
        rtask%p_rblockDiscrTrial => p_rblockDiscr
      end if
     
      ! Should we perform this task now?
      if (iand(rtask%iperform, iperform) .ne. 0) then
        ! Do we have block discretisations?
        if (.not.associated(rtask%p_rblockDiscrTest) .or.&
            .not.associated(rtask%p_rblockDiscrTrial)) then
          call output_line('Unable to find block discretisation',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
          call sys_halt()
        end if
        
        ! Clear preexisting block matrix
        call lsysbl_releaseMatrix(rmatrix)
        
        ! Create block matrix from discretisation
        if (associated(rtask%p_rblockDiscrTest,rtask%p_rblockDiscrTrial)) then
          call lsysbl_createMatBlockByDiscr(rtask%p_rblockDiscrTrial, rmatrix)
        else
          call lsysbl_createMatBlockByDiscr(rtask%p_rblockDiscrTrial,&
              rmatrix, rtask%p_rblockDiscrTest)
        end if
      end if

      !-------------------------------------------------------------------------
      ! Update submatrices
      !
      ! SYNTAX: submatrix(n) =
      !           matrixscalar@(i,j): MatrixSc1
      !           matrixscalar@(k,l): MatrixSc2
      !           ...
      !           matrixblock@(a,b):  MatrixBl1
      !
      ! For scalar matrices, the matrix objects MatrixSc1 and
      ! MatrixSc2 are generated at position (i,j) and (k,l) of the
      ! given block matrix, respectively. For block matrices, the
      ! (1,1) block of object MatrixBl1 is generated at position (a,b)
      ! of the given matrix.  The user is responsible to make sure
      ! that the given matrix provides all positions for the submatrices.

      call updateSubmatrix(rproblemLevel, rmatrix, rparlist,&
          ssectionName, 0, 0, iperform, p_rtaskList, rproblemTopLevel)
            
    elseif (trim(saction) .eq. 'DUPLICATE') then

      !-------------------------------------------------------------------------
      ! Duplicate block matrix or create 1-block wrapper for a scalar matrix
      !
      ! SYNTAX: matrixscalar = name,imatrix:#,ilev:#,...
      !     OR  matrixblock  = name,imatrix:#,ilev:#,...
      !
      ! If ilev is not given then the level of the current problem
      ! level structure is adopted.

      ! Find problem name from which to duplicate scalar matrix (if any)
      call parlst_getvalue_string(rparlist, ssectionName, 'matrixscalar', sparameter, '')
      call sys_toupper(sparameter)
      
      if (trim(sparameter) .ne. '') then
        
        !--- Scalar matrix case ------------------------------------------------

        ! Create new task for this scalar matrix
        nullify(rtask)
        allocate(rtask)
        nullify(rtask%p_rnextTask)
        
        ! Initialise optional parameters
        sproblem = trim(rproblemLevel%p_rproblem%cproblem)
        ilev     = rproblemLevel%ilev
        imatrix  = 1
        
        ! Get optional parameters if available
        call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
        
        do itoken = 1,ntoken
          call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
          
          iistart = 1
          call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
          call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
          
          if (trim(skeyword) .eq. 'PROBLEM') then
            sproblem = trim(svalue)
          elseif (trim(skeyword) .eq. 'ILEV') then
            read(svalue,'(I10)') ilev
          elseif (trim(skeyword) .eq. 'IMATRIX') then
            read(svalue,'(I10)') imatrix
          else
            call output_line('Unsupported keyword: '//trim(skeyword),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
            call sys_halt()
          end if
        end do
        
        ! Find problem level in global problem structure
        p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
        if (.not.associated(p_rproblemLevel)) then
          call output_line('Unable to find problem level in problem '//trim(sproblem),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
          call sys_halt()
        end if
        
        ! Get section name of scalar matrix
        call parlst_getvalue_string(rparlist, trim(sproblem),&
            'matrixscalar', ssubsectionName, isubstring=imatrix)
        
        ! Proceed with possibly prerequisite tasks
        call problem_updateMatrixScalar(p_rproblemLevel,&
            p_rproblemLevel%RmatrixScalar(imatrix), rparlist,&
            ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
        
        ! Specify new task for this block matrix
        rtask%bisMatrixScalar     = .false.
        rtask%ctask               =  PROBACTION_DUPLICATE
        rtask%ssectionName        =  ssectionName
        rtask%p_rmatrixScalarSrc  => p_rproblemLevel%RmatrixScalar(imatrix)
        rtask%p_rmatrixBlockDest  => rmatrix
        rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
        rtask%p_rproblemDest      => rproblemLevel%p_rproblem
        rtask%p_rproblemLevelSrc  => p_rproblemLevel
        rtask%p_rproblemLevelDest => rproblemLevel

      else

        !--- Block matrix case -------------------------------------------------

        ! Find problem name from which to duplicate entry of block matrix (if any)
        call parlst_getvalue_string(rparlist, ssectionName, 'matrixblock', sparameter, '')
        call sys_toupper(sparameter)
        
        if (trim(sparameter) .ne. '') then
          
          ! Create new task for this scalar matrix
          nullify(rtask)
          allocate(rtask)
          nullify(rtask%p_rnextTask)
          
          ! Initialise optional parameters
          sproblem = trim(rproblemLevel%p_rproblem%cproblem)
          ilev     = rproblemLevel%ilev
          imatrix  = 1

          ! Get optional parameters if available
          call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
          
          do itoken = 1,ntoken
            call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
            
            iistart = 1
            call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
            call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
            
            if (trim(skeyword) .eq. 'PROBLEM') then
              sproblem = trim(svalue)
            elseif (trim(skeyword) .eq. 'ILEV') then
              read(svalue,'(I10)') ilev
            elseif (trim(skeyword) .eq. 'IMATRIX') then
              read(svalue,'(I10)') imatrix
            else
              call output_line('Unsupported keyword: '//trim(skeyword),&
                  OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
              call sys_halt()
            end if
          end do

          ! Find problem level in global problem structure
          p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
          if (.not.associated(p_rproblemLevel)) then
            call output_line('Unable to find problem level in problem '//trim(sproblem),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
            call sys_halt()
          end if
          
          ! Get section name of block matrix
          call parlst_getvalue_string(rparlist, trim(sproblem),&
              'matrixblock', ssubsectionName, isubstring=imatrix)
          
          ! Proceed with possibly prerequisite tasks
          call problem_updateMatrixBlock(p_rproblemLevel,&
              p_rproblemLevel%RmatrixBlock(imatrix), rparlist,&
              ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)

          ! Specify new task for this block matrix
          rtask%bisMatrixScalar     = .false.
          rtask%ctask               =  PROBACTION_DUPLICATE
          rtask%ssectionName        =  ssectionName
          rtask%p_rmatrixBlockSrc   => p_rproblemLevel%RmatrixBlock(imatrix)
          rtask%p_rmatrixBlockDest  => rmatrix
          rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
          rtask%p_rproblemDest      => rproblemLevel%p_rproblem
          rtask%p_rproblemLevelSrc  => p_rproblemLevel
          rtask%p_rproblemLevelDest => rproblemLevel

        else
          call output_line('Either matrixscalar of matrixblock must be present',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
          call sys_halt()
        end if
      end if
    
      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this block matrix
          return
        end if
      else
        p_rtasklist => rtask
      end if      

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
          call sys_halt()
        end if
      end do

      ! Get optional parameters
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupstructure', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupStructure = get_dupType(sparameter)
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupcontent', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupContent = get_dupType(sparameter)

      ! Should we perform this task now?
      if (iand(rtask%iperform, iperform) .ne. 0) then
        if (associated(rtask%p_rmatrixScalarSrc)) then
          ! Duplicate scalar matrix as 1-block matrix
          call lsysbl_createMatFromScalar(rtask%p_rmatrixScalarSrc, rmatrixTemp)
          
          ! Copy content of temporal matrix by hand
          rmatrix%NEQ           = rmatrixTemp%NEQ
          rmatrix%NCOLS         = rmatrixTemp%NCOLS
          rmatrix%nblocksPerCol = 1
          rmatrix%nblocksPerRow = 1
          rmatrix%imatrixSpec   = LSYSBS_MSPEC_SCALAR

          allocate(rmatrix%RmatrixBlock(1,1))
          rmatrix%RmatrixBlock(1,1) = rmatrixTemp%RmatrixBlock(1,1)

          ! Release temporal matrix
          call lsysbl_releaseMatrix(rmatrixTemp)

        elseif(associated(rtask%p_rmatrixBlockSrc)) then
          ! Duplicate block matrix
          call lsysbl_duplicateMatrix(rtask%p_rmatrixBlockSrc,&
              rmatrix, rtask%cdupStructure, rtask%cdupContent)
        else
          call output_line('Neither scalar nor block matrix can be duplicated',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
          call sys_halt()
        end if
      end if
       
    else
      call output_line('Unsupported action: '//trim(saction),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
      call sys_halt()
    end if
    
    !---------------------------------------------------------------------------
    ! Second part: generate content
    !---------------------------------------------------------------------------

    nmethod = parlst_querysubstrings(rparlist, ssectionName, 'smethod')
    allocate(rtask%Cmethod(0:nmethod))

    ! What creation methods are adopted?
    do imethod = 0,nmethod

      call parlst_getvalue_string(rparlist, ssectionName,&
          'smethod', sparameter, '', isubstring=imethod)
      call sys_toupper(sparameter)
      
      if (trim(sparameter) .eq. 'VIRTUAL') then
        
        !-----------------------------------------------------------------------
        ! Create virtual matrix, aka, do nothing
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_VIRTUAL
        
      elseif (trim(sparameter) .eq. 'EMPTY') then
        
        !-----------------------------------------------------------------------
        ! Create empty matrix
        
        call lsysbl_allocEmptyMatrix(rmatrix, LSYSSC_SETM_ZERO)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_EMPTY
        
      elseif (trim(sparameter) .eq. 'IDENTITY') then
        
        !-----------------------------------------------------------------------
        ! Create identity matrix
        
        do j = 1,rmatrix%nblocksPerRow
          do i = 1,rmatrix%nblocksPerCol
            if (i.eq.j) then
              call lsyssc_initialiseIdentityMatrix(rmatrix%RmatrixBlock(i,j))
            else
              call lsyssc_allocEmptyMatrix(rmatrix%RmatrixBlock(i,j), LSYSSC_SETM_ZERO)
            end if
          end do
        end do
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_MATRIX_IDENTITY
        
      else
        call output_line('Unsupported method: '//trim(sparameter),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
        call sys_halt()
      end if
    end do

  contains

    ! Here, some auxiliary routines follow

    !***************************************************************************
    ! Search for given task in the list of tasks. If the task is not
    ! present in the list, then it is appended.
    function appendToTaskList(rtasklist, rtask) result(bexists)

      type(t_matrixTask), intent(inout), target :: rtasklist
      type(t_matrixTask), intent(in), target :: rtask
      
      ! local variable
      type(t_matrixTask), pointer :: p_rtask
      logical :: bexists      

      p_rtask => rtasklist
      do while(associated(p_rtask))
        if (associated(p_rtask%p_rmatrixBlockDest,&
            rtask%p_rmatrixBlockDest)) then
          bexists = .true.
          return
        end if
        if (.not.associated(p_rtask%p_rnextTask)) then
          p_rtask%p_rnextTask => rtask
          bexists = .false.
          return
        else
          p_rtask => p_rtask%p_rnextTask
        end if
      end do
      
    end function appendToTaskList

    !**************************************************************
    ! Return the numeric value of the duplication flad
    function get_dupType(sdupType) result(cdupType)

      character(len=*), intent(in) :: sdupType
      integer :: cdupType

      if (sdupType .eq. 'IGNORE') then
        cdupType = LSYSSC_DUP_IGNORE
      elseif (sdupType .eq. 'REMOVE') then
        cdupType = LSYSSC_DUP_REMOVE
      elseif (sdupType .eq. 'DISMISS') then
        cdupType = LSYSSC_DUP_DISMISS
      elseif (sdupType .eq. 'SHARE') then
        cdupType = LSYSSC_DUP_SHARE
      elseif (sdupType .eq. 'COPY') then
        cdupType = LSYSSC_DUP_COPY
      elseif (sdupType .eq. 'COPYOVERWRITE') then
        cdupType = LSYSSC_DUP_COPYOVERWRITE
      elseif (sdupType .eq. 'ASIS') then
        cdupType = LSYSSC_DUP_ASIS
      elseif (sdupType .eq. 'EMPTY') then
        cdupType = LSYSSC_DUP_EMPTY
      elseif (sdupType .eq. 'TEMPLATE') then
        cdupType = LSYSSC_DUP_TEMPLATE
      else
        cdupType = 0
      end if

    end function get_dupType
    
    !**************************************************************
    ! Update the submatrix
    subroutine updateSubmatrix(rproblemLevel, rmatrix,&
        rparlist, ssectionName, iblockrowOffset, iblockcolOffset,&
        iperform, p_rtaskList, rproblemTopLevel)

      type(t_problemLevel), intent(in), target :: rproblemLevel
      type(t_problem), intent(in), target, optional :: rproblemTopLevel
      type(t_matrixBlock), intent(inout) :: rmatrix
      type(t_parlist), intent(in) :: rparlist
      character(len=*), intent(in) :: ssectionName
      integer, intent(in) :: iblockrowOffset,iblockcolOffset
      integer(i32), intent(in) :: iperform
      type(t_matrixTask), pointer :: p_rtaskList

      ! local variables
      type(t_matrixBlock) :: rmatrixTemp
      character(len=PARLST_MLDATA) :: sparameter,stoken,ssubtoken,ssubmatrix
      integer :: i,iblockcol,iblockrow,iiistart,iistart,istart,isubmatrix,&
                 j,nsubmatrix

      ! Get number of submatrices
      nsubmatrix = parlst_querysubstrings(rparlist, ssectionName, 'submatrix')

      ! Loop over all submatrices
      do isubmatrix=1,nsubmatrix
        
        call parlst_getvalue_string(rparlist, ssectionName, 'submatrix',&
            sparameter, isubstring=isubmatrix)
        call sys_toupper(sparameter)

        ! Extract content from parameter file line by line
        ! SYNTAX: matrixscalar@(1,1): MatrixSc1
        !     OR  matrixblock@(1,1):  MatrixBl1
        istart = 1; iistart = 1; iiistart = 1
        call sys_getNextToken (sparameter, stoken, istart, ":", .false.)
        call sys_getNextToken (stoken, ssubmatrix, iistart, "@", .false.)
        call sys_getNextToken (stoken, ssubtoken, iiistart, ",", .false.)
        read(ssubtoken(iistart:iiistart),*) iblockrow
        iistart = iiistart
        call sys_getNextToken (stoken, ssubtoken, iiistart, ")", .false.)
        read(ssubtoken,*) iblockcol
        call sys_getNextToken (sparameter, stoken, istart, ":", .false.)

        ! Apply given block offsets
        iblockcol = iblockcol+iblockcolOffset
        iblockrow = iblockrow+iblockrowOffset

        ! Check of this is a scalar matrix or a block matrix
        if (trim(ssubmatrix) .eq. "MATRIXSCALAR") then
          
          ! Check if scalar submatrix is available
          if ((rmatrix%nblocksPerCol .lt. iblockrow) .or.&
              (rmatrix%nblocksPerRow .lt. iblockcol)) then
            call output_line('Position of submatrix is not available',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
            call sys_halt()
          end if

          ! Update single scalar submatrix
          call problem_updateMatrixScalar(rproblemLevel,&
              rmatrix%RmatrixBlock(iblockrow,iblockcol), rparlist,&
              trim(stoken), p_rtasklist, rproblemTopLevel, iperform)
          
        elseif(trim(ssubmatrix) .eq. "MATRIXBLOCK") then
          
          ! Update single block submatrix
          call problem_updateMatrixBlock(rproblemLevel,&
              rmatrixTemp, rparlist, trim(stoken), p_rtasklist,&
              rproblemTopLevel, iperform)

          ! Transfer content of scalar submatrices
          do j=1,rmatrixTemp%nblocksPerRow
            do i=1,rmatrixTemp%nblocksPerCol
              rmatrix%RmatrixBlock(iblockrow+i-1,iblockcol+j-1)=&
                  rmatrixTemp%RmatrixBlock(i,j)
            end do
          end do

          ! Release temporal matrix
          call lsysbl_releaseMatrix(rmatrixTemp)
          
        else
          call output_line('Either matrixscalar of matrixblock must be present',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrixBlock')
          call sys_halt()
        end if
        
      end do

    end subroutine updateSubmatrix

  end subroutine problem_updateMatrixBlock

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_updateMatrix(rtasklist, iperformSpec)

!<description>
    ! This subroutine updates the scalar/block matrices with the
    ! internal values assigned to the items of the task list.
!</description>

!<input>
    ! task list
    type(t_matrixTask), intent(in), target :: rtasklist

    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>
!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type(t_matrixBlock) :: rmatrixTemp
    type(t_matrixTask), pointer :: p_rtask
    integer :: iperform,imethod

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec
    
    ! Iterate over all tasks
    p_rtask => rtasklist
    taskloop: do while (associated(p_rtask))

      ! Do we have to perform this task?
      if (iand(p_rtask%iperform, iperform) .eq. 0) then
        p_rtask => p_rtask%p_rnextTask
        cycle taskloop
      end if
      
      !-------------------------------------------------------------------------
      ! First part: create or duplicate matrix
      !-------------------------------------------------------------------------
      select case(p_rtask%ctask)
      case(PROBACTION_CREATE)

        ! Do we have to assemble a scalar matrix?
        if (p_rtask%bisMatrixScalar) then
          
          !---------------------------------------------------------------------
          ! Create scalar matrix
          
          ! Do we have spatial discretisations?
          if (.not.associated(p_rtask%p_rspatialDiscrTest) .or.&
              .not.associated(p_rtask%p_rspatialDiscrTrial)) then
            call output_line('Unable to find spatial discretisation',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrix')
            call sys_halt()
          end if
          
          ! Clear preexisting matrix
          call lsyssc_releaseMatrix(p_rtask%p_rmatrixScalarDest)

          ! Create matrix structure and set data type
          if (associated(p_rtask%p_rspatialDiscrTest,&
                         p_rtask%p_rspatialDiscrTrial)) then
            call bilf_createMatrixStructure(p_rtask%p_rspatialDiscrTrial,&
                p_rtask%cmatrixFormat, p_rtask%p_rmatrixScalarDest)
          else
            call bilf_createMatrixStructure(p_rtask%p_rspatialDiscrTrial,&
                p_rtask%cmatrixFormat, p_rtask%p_rmatrixScalarDest,&
                p_rtask%p_rspatialDiscrTest)
          end if
          call lsyssc_setDataTypeMatrix (p_rtask%p_rmatrixScalarDest, p_rtask%cdataType)
          p_rtask%p_rmatrixScalarDest%dscaleFactor = p_rtask%dscaleFactor
          if ((p_rtask%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) .or.&
              (p_rtask%cmatrixFormat .eq. LSYSSC_MATRIX9INTL)) then
            p_rtask%p_rmatrixScalarDest%nvar = p_rtask%nvar
          end if

        else

          !---------------------------------------------------------------------
          ! Create block matrix

          ! Do we have block discretisations?
          if (.not.associated(p_rtask%p_rblockDiscrTest) .or.&
              .not.associated(p_rtask%p_rblockDiscrTrial)) then
            call output_line('Unable to find block discretisation',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrix')
            call sys_halt()
          end if

          ! Clear preexisting block matrix
          call lsysbl_releaseMatrix(p_rtask%p_rmatrixBlockDest)

          ! Create block matrix from discretisation
          if (associated(p_rtask%p_rblockDiscrTest,p_rtask%p_rblockDiscrTrial)) then
            call lsysbl_createMatBlockByDiscr(p_rtask%p_rblockDiscrTrial,&
                p_rtask%p_rmatrixBlockDest)
          else
            call lsysbl_createMatBlockByDiscr(p_rtask%p_rblockDiscrTrial,&
                p_rtask%p_rmatrixBlockDest, p_rtask%p_rblockDiscrTest)
          end if
          
          ! Proceed with next task
          p_rtask => p_rtask%p_rnextTask
          
          ! That`s it for block matrices
          cycle taskloop
        end if

      case(PROBACTION_DUPLICATE)

        !-----------------------------------------------------------------------
        ! Duplicate scalar/block matrix

        if (p_rtask%bisMatrixScalar) then
          ! Duplicate scalar matrix
          call lsyssc_duplicateMatrix(p_rtask%p_rmatrixScalarSrc,&
              p_rtask%p_rmatrixScalarDest, p_rtask%cdupStructure, p_rtask%cdupContent)
        else
          if (associated(p_rtask%p_rmatrixScalarSrc)) then
            ! Duplicate scalar matrix as 1-block matrix
            call lsysbl_createMatFromScalar(p_rtask%p_rmatrixScalarSrc, rmatrixTemp)

            ! Copy content of temporal matrix by hand
            p_rtask%p_rmatrixBlockDest%NEQ           = rmatrixTemp%NEQ
            p_rtask%p_rmatrixBlockDest%NCOLS         = rmatrixTemp%NCOLS
            p_rtask%p_rmatrixBlockDest%nblocksPerCol = 1
            p_rtask%p_rmatrixBlockDest%nblocksPerRow = 1
            p_rtask%p_rmatrixBlockDest%imatrixSpec   = LSYSBS_MSPEC_SCALAR

            allocate(p_rtask%p_rmatrixBlockDest%RmatrixBlock(1,1))
            p_rtask%p_rmatrixBlockDest%RmatrixBlock(1,1) = rmatrixTemp%RmatrixBlock(1,1)

            ! Release temporal matrix
            call lsysbl_releaseMatrix(rmatrixTemp)

          elseif(associated(p_rtask%p_rmatrixBlockSrc)) then
            ! Duplicate block matrix
            call lsysbl_duplicateMatrix(p_rtask%p_rmatrixBlockSrc,&
                p_rtask%p_rmatrixBlockDest, p_rtask%cdupStructure, p_rtask%cdupContent)
          else
            call output_line('Neither scalar nor block matrix can be duplicated',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrix')
            call sys_halt()
          end if
        end if

      case default
        call output_line('Unsupported action: '//sys_siL(p_rtask%ctask,3),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMamtrix')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------
      ! Second part: generate content (only for scalar matrices)
      !-------------------------------------------------------------------------

      ! What creation methods are adopted?
      do imethod = lbound(p_rtask%Cmethod,1), ubound(p_rtask%Cmethod,1)

        select case (p_rtask%Cmethod(imethod))
        case(PROBMETH_MATRIX_VIRTUAL)
          
          !---------------------------------------------------------------------
          ! Create virtual matrix, aka, do nothing
          
        case(PROBMETH_MATRIX_EMPTY)
          
          !---------------------------------------------------------------------
          ! Create empty matrix
          
          call lsyssc_allocEmptyMatrix(p_rtask%p_rmatrixScalarDest, LSYSSC_SETM_ZERO)
          
        case(PROBMETH_MATRIX_IDENTITY)
          
          !---------------------------------------------------------------------
          ! Create identity matrix
          
          call lsyssc_initialiseIdentityMatrix(p_rtask%p_rmatrixScalarDest)
          
        case(PROBMETH_MATRIX_BILF)
          
          !---------------------------------------------------------------------
          ! Create matrix by bilinearform evaluation         
          
          if (p_rtask%p_rform%ballCoeffConstant) then
            if (associated(p_rtask%p_rcubatureInfo)) then
              call bilf_buildMatrixScalar(p_rtask%p_rform, .true.,&
                  p_rtask%p_rmatrixScalarDest, p_rtask%p_rcubatureInfo)
            else
              call bilf_buildMatrixScalar(p_rtask%p_rform, .true.,&
                  p_rtask%p_rmatrixScalarDest)
            end if
          else
            rcollection%p_rfparserQuickAccess1 => p_rtask%p_rfparser
            if (associated(p_rtask%p_rcubatureInfo)) then
              call bilf_buildMatrixScalar(p_rtask%p_rform, .true.,&
                  p_rtask%p_rmatrixScalarDest, p_rtask%p_rcubatureInfo,&
                  problem_buildMatrixSc_fparser_sim, rcollection)
            else
              call bilf_buildMatrixScalar(p_rtask%p_rform, .true.,&
                  p_rtask%p_rmatrixScalarDest,&
                  problem_buildMatrixSc_fparser_sim, rcollection)
            end if
          end if
          
        case(PROBMETH_MATRIX_LUMP_STD)
          
          !---------------------------------------------------------------------
          ! Create matrix by standard lumping
          
          call lsyssc_lumpMatrixScalar(p_rtask%p_rmatrixScalarDest, LSYSSC_LUMP_STD)
          
        case(PROBMETH_MATRIX_LUMP_DIAG)
          
          !---------------------------------------------------------------------
          ! Create matrix by diagonal lumping
          
          call lsyssc_lumpMatrixScalar(p_rtask%p_rmatrixScalarDest, LSYSSC_LUMP_DIAG)
          
        case default
          call output_line('Unsupported method: '//sys_siL(p_rtask%Cmethod(imethod),3),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMatrix')
          call sys_halt()
        end select
      end do

      ! Proceed with next task
      p_rtask => p_rtask%p_rnextTask
    end do taskloop
    
  end subroutine problem_updateMatrix

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_clearMatrixTL(p_rtasklist)

!<description>
    ! This subroutine clears the scalar/block matrix task list.
!</description>

!<input>
    ! task list
    type(t_matrixTask), pointer :: p_rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_matrixTask), pointer :: p_rtask,p_rtaskSave

    p_rtask => p_rtasklist
    do while(associated(p_rtask))
      if (associated(p_rtask%Cmethod)) deallocate(p_rtask%Cmethod)
      if (associated(p_rtask%p_rform)) deallocate(p_rtask%p_rform)
      if (associated(p_rtask%p_rfparser)) then
        call fparser_release(p_rtask%p_rfparser)
        deallocate(p_rtask%p_rfparser)
      end if
      p_rtaskSave => p_rtask
      p_rtask     => p_rtask%p_rnextTask
      deallocate(p_rtaskSave)
    end do
    
  end subroutine problem_clearMatrixTL

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_infoMatrixTL(rtasklist)

!<description>
    ! This subroutine outputs information about the scalar/block
    ! matrix task list.
!</description>

!<input>
    ! task list
    type(t_matrixTask), intent(in), target :: rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_matrixTask), pointer :: p_rtask
    integer :: i

    p_rtask => rtasklist
    do while(associated(p_rtask))

      select case(p_rtask%ctask)
      case (PROBACTION_NONE)
        call output_lbrk()
        call output_line ('MatrixTask: DO NOTHING')

      case (PROBACTION_CREATE)
        call output_lbrk()
        call output_line ('MatrixTask: CREATE')

      case (PROBACTION_DUPLICATE)
        call output_lbrk()
        call output_line ('MatrixTask: DUPLICATE')

      case default
        call output_line ('MatrixTask: UNSUPPORTED')
      end select

      call output_line ('-----------------------')
      call output_line ('matrix type               : '//merge('SCALAR', 'BLOCK ', p_rtask%bisMatrixScalar))
      call output_line ('section name              : '//trim(p_rtask%ssectionName))
      call output_line ('destination problem       : '//trim(p_rtask%p_rproblemDest%cproblem))
      call output_line ('destination problem level : '//trim(sys_siL(p_rtask%p_rproblemLevelDest%ilev,15)))
      call output_line ('destination matrix scalar : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rmatrixScalarDest)))
      call output_line ('destination matrix block  : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rmatrixBlockDest)))
      if (associated(p_rtask%p_rproblemSrc))&
      call output_line ('source problem            : '//trim(p_rtask%p_rproblemSrc%cproblem))
      if (associated(p_rtask%p_rproblemLevelSrc))&
      call output_line ('source problem level      : '//trim(sys_siL(p_rtask%p_rproblemLevelSrc%ilev,15)))
      call output_line ('source matrix scalar      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rmatrixScalarSrc)))
      call output_line ('source matrix block       : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rmatrixBlockSrc)))
      call output_line ('iperform                  : '//trim(sys_siL(int(p_rtask%iperform),15)))
      call output_line ('method                    : '//trim(sys_siL(p_rtask%Cmethod(0),15)))
      do i = 1,ubound(p_rtask%Cmethod,1)
        call output_line ('                          : '//trim(sys_siL(p_rtask%Cmethod(i),15)))
      end do
      call output_line ('data type                 : '//trim(sys_siL(p_rtask%cdatatype,15)))
      call output_line ('matrix format             : '//trim(sys_siL(p_rtask%cmatrixFormat,15)))
      call output_line ('matrix interleave format  : '//trim(sys_siL(p_rtask%cinterleavematrixFormat,15)))
      call output_line ('nvar                      : '//trim(sys_siL(p_rtask%nvar,15)))
      call output_line ('scaling factor            : '//trim(sys_sdL(p_rtask%dscaleFactor,15)))
      call output_line ('spatial discr. trial space: '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rspatialDiscrTrial)))
      call output_line ('spatial discr. test space : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rspatialDiscrTest)))
      call output_line ('block discr. trial space  : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rblockDiscrTrial)))
      call output_line ('block discr. test space   : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rblockDiscrTest)))
      call output_line ('cubature info structure   : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rcubatureInfo)))
      call output_line ('function parser           : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rfparser)))
      call output_line ('bilinear form evaluator   : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rform)))

      p_rtask => p_rtask%p_rnextTask
    end do
    
  end subroutine problem_infoMatrixTL
  
  !*****************************************************************************

!<subroutine>

  subroutine problem_initVectorAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises all scalar/block vectors of the
    ! given problem level with the values provided by the parameter
    ! list.  If the optional parameter p_rtasklist is given, then the
    ! task list which is created during the initialisation is returned.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! OPTIONAL: task list
    type(t_vectorTask), intent(inout), pointer, optional :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_vectorTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update all scalar/block vectors
    call problem_updateVectorAll(rproblemLevel,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearVectorTL(p_rtasklistLocal)
    end if

  end subroutine problem_initVectorAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_updateVectorAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine the task list for all scalar/block vectors of the
    ! given problem level with the values provided by the parameter list.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! pointer to the task list
    type(t_vectorTask), pointer :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: ssubsectionName
    integer :: i,nvectorScalar,nvectorBlock

    ! Consistency check
    nvectorScalar = parlst_querysubstrings(rparlist, ssectionName,&
                                           'vectorscalar')
    if ((nvectorScalar .gt. 0) .and.&
        (nvectorScalar .ne. size(rproblemLevel%RvectorScalar))) then
      call output_line('Invalid number of scalar vectors',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateAllVectorTaskListFromParlst')
      call sys_halt()
    end if

    nvectorBlock = parlst_querysubstrings(rparlist, ssectionName,&
                                          'vectorblock')
    if ((nvectorBlock .gt. 0) .and.&
        (nvectorBlock .ne. size(rproblemLevel%RvectorBlock))) then
      call output_line('Invalid number of block vectors',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateAllVectorTaskListFromParlst')
      call sys_halt()
    end if
    
    ! Update all scalar vectors
    do i=1,nvectorScalar
      call parlst_getvalue_string(rparlist, ssectionName,&
          'vectorscalar', ssubsectionName, isubstring=i)
      call problem_updateVectorScalar(rproblemLevel,&
          rproblemLevel%RvectorScalar(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do
    
    ! Update all block vectors
    do i=1,nvectorBlock
      call parlst_getvalue_string(rparlist, ssectionName,&
          'vectorblock', ssubsectionName, isubstring=i)
      call problem_updateVecBl(rproblemLevel,&
          rproblemLevel%RvectorBlock(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do

  end subroutine problem_updateVectorAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_initVectorScalar(rproblemLevel, rvector,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given scalar vector on the given
    ! problem level with the values provided by the parameter list. If
    ! the optional parameter p_rtasklist is given, then the task list
    ! which is created during the initialisation is returned.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar/block
    ! vectors need to be created and copied
    type(t_vectorTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! scalar vector
    type(t_vectorScalar), intent(out), target :: rvector
!</output>
!</subroutine>

    ! local variables
    type(t_vectorTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given scalar vector
    call problem_updateVectorScalar(rproblemLevel, rvector,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearVectorTL(p_rtasklistLocal)
    end if

  end subroutine problem_initVectorScalar

  !*****************************************************************************

!<subroutine>

  subroutine problem_initVectorBlock(rproblemLevel, rvector,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given block vector on the given
    ! problem level with the values provided by the parameter list. If
    ! the optional parameter p_rtasklist is given, then the task list
    ! which is created during the initialisation is returned.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar/block
    ! vectors need to be created and copied
    type(t_vectorTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! block vector
    type(t_vectorBlock), intent(out), target :: rvector
!</output>
!</subroutine>

    ! local variables
    type(t_vectorTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given block vector
    call problem_updateVecBl(rproblemLevel, rvector,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearVectorTL(p_rtasklistLocal)
    end if

  end subroutine problem_initVectorBlock

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateVectorScalar(rproblemLevel,&
      rvector, rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given scalar vector on the given problem level with the values
    ! provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar vectors
    ! need to be created and copied
    type(t_vectorTask), pointer :: p_rtasklist

    ! scalar vector
    type(t_vectorScalar), intent(inout), target :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvectorBlock
    type(t_vectorTask), pointer :: rtask
    type(t_problem), pointer :: p_rproblemTopLevel
    type(t_problemLevel), pointer :: p_rproblemLevel
    character(len=PARLST_MLDATA) :: skeyword,sparameter,sproblem,stoken,svalue
    character(len=SYS_STRLEN) :: ssubsectionName
    character(len=SYS_STRLEN) :: saction,sperform,sperf
    integer :: iblock,iistart,ilev,ivector,istart,itoken,ntoken,imethod,nmethod
    integer(i32) :: iperform

    ! Set pointer to problem structure
    p_rproblemTopLevel => rproblemLevel%p_rproblem
    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec

    ! Get type of action
    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
    call sys_toupper(saction)
    call sys_toupper(sperform)

    !---------------------------------------------------------------------------
    ! First part: generate structure
    !---------------------------------------------------------------------------

    ! What action should be performed?
    if (trim(saction) .eq. 'NONE') then
      !-------------------------------------------------------------------------
      ! Do nothing
      
    elseif (trim(saction) .eq. 'CREATE') then
      
      !-------------------------------------------------------------------------
      ! Create scalar vector
      !
      ! SYNTAX: vectorscalar(n) =
      !           VectorScalar1
      !                ...
      !           VectorScalar2

      ! Create new task
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%bisVectorScalar     = .true.
      rtask%ctask               =  PROBACTION_CREATE
      rtask%ssectionName        =  ssectionName
      rtask%p_rvectorScalarDest => rvector
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelDest => rproblemLevel

      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this scalar vector
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
          call sys_halt()
        end if
      end do
      
      ! Get optional parameters
      call parlst_getvalue_int(rparlist, ssectionName,&
          'cdatatype', rtask%cdataType, ST_NOHANDLE)
      if (rtask%cdataType .eq. ST_NOHANDLE) then
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sdatatype', sparameter, 'DOUBLE')
        call sys_toupper(sparameter)
        rtask%cdataType = get_dataType(trim(sparameter))
      end if
      
      call parlst_getvalue_int(rparlist, ssectionName,&
          'nvar', rtask%nvar, 1)
      
      ! Get spatial discretisation for test and trial space
      call parlst_getvalue_string(rparlist, ssectionName,&
          'discretisation', sparameter)
      call sys_toupper(sparameter)
      rtask%p_rspatialDiscr => problem_getSpatialDiscr(&
          rproblemLevel%p_rproblem, rproblemLevel, sparameter)
      
      ! Should we perform this task now?
      if (iand(rtask%iperform, iperform) .ne. 0) then
        ! Do we have spatial discretisations?
        if (.not.associated(rtask%p_rspatialDiscr)) then
          call output_line('Unable to find spatial discretisation',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
          call sys_halt()
        end if
        
        ! Clear preexisting scalar vector
        call lsyssc_releaseVector(rvector)
        
        ! Create vector structure and set data type
        call lsyssc_createVector(rtask%p_rspatialDiscr, rvector,&
            rtask%nvar, .false., rtask%cdataType)
      end if

    elseif (trim(saction) .eq. 'DUPLICATE') then

      !-------------------------------------------------------------------------
      ! Duplicate scalar vector or subvector of a block vector
      !
      ! SYNTAX: vectorscalar = name,ivector:#,ilev:#,...
      !     OR  vectorblock  = name,ivector:#,ilev:#,iblock:#
      !
      ! If ilev is not given then the level of the current problem
      ! level structure is adopted. If iblock are not given then
      ! standard value 1 is used in both cases

      ! Find problem name from which to duplicate scalar vector (if any)
      call parlst_getvalue_string(rparlist, ssectionName, 'vectorscalar', sparameter, '')
      call sys_toupper(sparameter)

      if (trim(sparameter) .ne. '') then
        
        !--- Scalar vector case ------------------------------------------------

        ! Create new task for this scalar vector
        nullify(rtask)
        allocate(rtask)
        nullify(rtask%p_rnextTask)
        
        ! Initialise optional parameters
        sproblem = trim(rproblemLevel%p_rproblem%cproblem)
        ilev     = rproblemLevel%ilev
        ivector  = 1
        
        ! Get optional parameters if available
        call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
        
        do itoken = 1,ntoken
          call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
          
          iistart = 1
          call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
          call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
          
          if (trim(skeyword) .eq. 'PROBLEM') then
            sproblem = trim(svalue)
          elseif (trim(skeyword) .eq. 'ILEV') then
            read(svalue,'(I10)') ilev
          elseif (trim(skeyword) .eq. 'IVECTOR') then
            read(svalue,'(I10)') ivector
          else
            call output_line('Unsupported keyword: '//trim(skeyword),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
            call sys_halt()
          end if
        end do
        
        ! Find problem level in global problem structure
        p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
        if (.not.associated(p_rproblemLevel)) then
          call output_line('Unable to find problem level in problem '//trim(sproblem),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
          call sys_halt()
        end if
        
        ! Get section name of scalar vector
        call parlst_getvalue_string(rparlist, trim(sproblem),&
            'vectorscalar', ssubsectionName, isubstring=ivector)

        ! Proceed with possibly prerequisite tasks
        call problem_updateVectorScalar(p_rproblemLevel,&
            p_rproblemLevel%RvectorScalar(ivector), rparlist,&
            ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
        
        ! Specify new task for this scalar vector
        rtask%bisVectorScalar     = .true.
        rtask%ctask               =  PROBACTION_DUPLICATE
        rtask%ssectionName        =  ssectionName
        rtask%p_rvectorScalarSrc  => p_rproblemLevel%RvectorScalar(ivector)
        rtask%p_rvectorScalarDest => rvector
        rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
        rtask%p_rproblemDest      => rproblemLevel%p_rproblem
        rtask%p_rproblemLevelSrc  => p_rproblemLevel
        rtask%p_rproblemLevelDest => rproblemLevel

      else

        !--- Block vector case -------------------------------------------------

        ! Find problem name from which to duplicate entry of block vector (if any)
        call parlst_getvalue_string(rparlist, ssectionName, 'vectorblock', sparameter, '')
        call sys_toupper(sparameter)
      
        if (trim(sparameter) .ne. '') then
          
          ! Create new task for this scalar vector
          nullify(rtask)
          allocate(rtask)
          nullify(rtask%p_rnextTask)
          
          ! Initialise optional parameters
          sproblem = trim(rproblemLevel%p_rproblem%cproblem)
          ilev     = rproblemLevel%ilev
          ivector  = 1
          iblock   = 1
          
          ! Get optional parameters if available
          call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
          
          do itoken = 1,ntoken
            call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
            
            iistart = 1
            call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
            call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
            
            if (trim(skeyword) .eq. 'PROBLEM') then
              sproblem = trim(svalue)
            elseif (trim(skeyword) .eq. 'ILEV') then
              read(svalue,'(I10)') ilev
            elseif (trim(skeyword) .eq. 'IVECTOR') then
              read(svalue,'(I10)') ivector
            elseif (trim(skeyword) .eq. 'IBLOCK') then
              read(svalue,'(I10)') iblock
            else
              call output_line('Unsupported keyword: '//trim(skeyword),&
                  OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
              call sys_halt()
            end if
          end do
          
          ! Find problem level in global problem structure
          p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
          if (.not.associated(p_rproblemLevel)) then
            call output_line('Unable to find problem level in problem '//trim(sproblem),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
            call sys_halt()
          end if
          
          ! Get section name of block vector
          call parlst_getvalue_string(rparlist, trim(sproblem),&
              'vectorblock', ssubsectionName, isubstring=ivector)

          ! Proceed with possibly prerequisite tasks
          call problem_updateVecBl(p_rproblemLevel,&
              p_rproblemLevel%RvectorBlock(ivector), rparlist,&
              ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)

          ! Check if scalar subvector is available
          p_rvectorBlock => p_rproblemLevel%RvectorBlock(ivector)
          if ((size(p_rvectorBlock%RvectorBlock,1) .lt. iblock)) then
            call output_line('Scalar subvector of block vector is not available',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
            call sys_halt()
          end if

          ! Specify new task for this scalar vector
          rtask%bisVectorScalar     = .true.
          rtask%ctask               =  PROBACTION_DUPLICATE
          rtask%ssectionName        =  ssectionName
          rtask%p_rvectorScalarSrc  => p_rvectorBlock%RvectorBlock(iblock)
          rtask%p_rvectorScalarDest => rvector
          rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
          rtask%p_rproblemDest      => rproblemLevel%p_rproblem
          rtask%p_rproblemLevelSrc  => p_rproblemLevel
          rtask%p_rproblemLevelDest => rproblemLevel

        else
          call output_line('Either vectorscalar of vectorblock must be present',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
          call sys_halt()
        end if
      end if

      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this scalar vector
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
          call sys_halt()
        end if
      end do
          
      ! Get optional parameters
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupstructure', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupStructure = get_dupType(sparameter)
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupcontent', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupContent = get_dupType(sparameter)

      ! Duplicate scalar vector
      if (iand(rtask%iperform, iperform) .ne. 0) then
        call lsyssc_duplicateVector(rtask%p_rvectorScalarSrc,&
            rvector, rtask%cdupStructure, rtask%cdupContent)
      end if
      
    else
      call output_line('Unsupported action: '//trim(saction),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
      call sys_halt()
    end if
    
    !---------------------------------------------------------------------------
    ! Second part: generate content
    !---------------------------------------------------------------------------

    nmethod = parlst_querysubstrings(rparlist, ssectionName, 'smethod')
    allocate(rtask%Cmethod(0:nmethod))

    ! What creation methods are adopted?
    do imethod = 0,nmethod
      
      call parlst_getvalue_string(rparlist, ssectionName,&
          'smethod', sparameter, '', isubstring=imethod)
      call sys_toupper(sparameter)

      
      if (trim(sparameter) .eq. 'VIRTUAL') then
        
        !-----------------------------------------------------------------------
        ! Create virtual vector, aka, do nothing
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_VIRTUAL

      elseif (trim(sparameter) .eq. 'EMPTY') then
        
        !-----------------------------------------------------------------------
        ! Create empty vector
        
        call lsyssc_clearVector(rvector, 0.0_DP)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_EMPTY
      
      elseif (trim(sparameter) .eq. 'UNITY') then
        
        !-----------------------------------------------------------------------
        ! Create unit vector
        
        call lsyssc_clearVector(rvector, 1.0_DP)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_UNITY
      
      elseif (trim(sparameter) .eq. 'LINF') then
        
        !-----------------------------------------------------------------------
        ! Create vector by linearform evaluation         
        
        call assembleScalarVectorLinf(rparlist, rtask)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_LINF
      
      elseif (trim(sparameter) .eq. 'DOF') then
        
        !-----------------------------------------------------------------------
        ! Create coordinate vector for degrees of freedom
        
        call lin_calcDofCoords(rtask%p_rspatialDiscr, rvector)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_DOF
        
      else
        call output_line('Unsupported method: '//trim(sparameter),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
        call sys_halt()
      end if
    end do

  contains

    ! Here, some auxiliary routines follow

    !***************************************************************************
    ! Search for given task in the list of tasks. If the task is not
    ! present in the list, then it is appended.
    function appendToTaskList(rtasklist, rtask) result(bexists)

      type(t_vectorTask), intent(inout), target :: rtasklist
      type(t_vectorTask), intent(in), target :: rtask
      
      ! local variable
      type(t_vectorTask), pointer :: p_rtask
      logical :: bexists      

      p_rtask => rtasklist
      do while(associated(p_rtask))
        if (associated(p_rtask%p_rvectorScalarDest,&
            rtask%p_rvectorScalarDest)) then
          bexists = .true.
          return
        end if
        if (.not.associated(p_rtask%p_rnextTask)) then
          p_rtask%p_rnextTask => rtask
          bexists = .false.
          return
        else
          p_rtask => p_rtask%p_rnextTask
        end if
      end do
      
    end function appendToTaskList

    !**************************************************************
    ! Return the numeric value of the data type
    function get_dataType(sdataType) result(cdataType)

      character(len=*), intent(in) :: sdataType
      integer :: cdataType

      if (sdataType .eq. 'QUAD') then
        cdataType = ST_QUAD
      elseif (sdataType .eq. 'DOUBLE') then
        cdataType = ST_DOUBLE
      elseif (sdataType .eq. 'SINGLE') then
        cdataType = ST_SINGLE
      else
        cdataType = ST_NOHANDLE
      end if

    end function get_dataType

    !**************************************************************
    ! Return the numeric value of the duplication flad
    function get_dupType(sdupType) result(cdupType)

      character(len=*), intent(in) :: sdupType
      integer :: cdupType

      if (sdupType .eq. 'IGNORE') then
        cdupType = LSYSSC_DUP_IGNORE
      elseif (sdupType .eq. 'REMOVE') then
        cdupType = LSYSSC_DUP_REMOVE
      elseif (sdupType .eq. 'DISMISS') then
        cdupType = LSYSSC_DUP_DISMISS
      elseif (sdupType .eq. 'SHARE') then
        cdupType = LSYSSC_DUP_SHARE
      elseif (sdupType .eq. 'COPY') then
        cdupType = LSYSSC_DUP_COPY
      elseif (sdupType .eq. 'COPYOVERWRITE') then
        cdupType = LSYSSC_DUP_COPYOVERWRITE
      elseif (sdupType .eq. 'ASIS') then
        cdupType = LSYSSC_DUP_ASIS
      elseif (sdupType .eq. 'EMPTY') then
        cdupType = LSYSSC_DUP_EMPTY
      elseif (sdupType .eq. 'TEMPLATE') then
        cdupType = LSYSSC_DUP_TEMPLATE
      else
        cdupType = 0
      end if

    end function get_dupType
    
    !**************************************************************
    ! Assemble the scalar vector from linear form evaluation
    subroutine assembleScalarVectorLinf(rparlist, rtask)

      type(t_parlist), intent(in) :: rparlist
      type(t_vectorTask), intent(inout) :: rtask

      ! local variables
      type(t_collection) :: rcollection
      type(t_scalarCubatureInfo), pointer :: p_rcubInfo
      character(len=PARLST_MLDATA) :: sparameter,stoken
      integer :: i,j,itermCount,istart

      ! Get cubature info structure
      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
          'cubatureinfo', sparameter)
      call sys_toupper(sparameter)
      p_rcubInfo => problem_getCubInfo(rtask%p_rproblemDest,&
          rtask%p_rproblemLevelDest, sparameter)
      
      ! Update internal data of the task item
      rtask%p_rcubatureInfo => p_rcubInfo
      allocate(rtask%p_rform)
      
      ! Setup trial functions
      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
          'sfunction', sparameter)
      call sys_countTokens(sparameter, itermCount, ',', .false.)
      rtask%p_rform%itermCount = itermCount
      
      istart = 1
      do i = 1,itermCount
        call sys_getNextToken(sparameter, stoken, istart, ',', .false.)
        rtask%p_rform%Idescriptors(i) = der_igetID(stoken)
      end do

      ! Setup coefficients
      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
          'scoefficient', sparameter, '')
      call sys_countTokens(sparameter, itermCount, ',', .false.)
      
      if (itermCount .gt. 0) then
        ! Use individual coefficients for each component
        rtask%p_rform%itermCount = min(rtask%p_rform%itermCount,itermCount)
        
        istart = 1
        do i=1,itermCount
          call sys_getNextToken(sparameter, stoken, istart, ',', .false.)
          if (sys_isNumeric(stoken)) then
            read(stoken,*) rtask%p_rform%Dcoefficients(i)
          else
            ! Not all coefficients are constant. Thus, we create a
            ! function parser object and perform bilinear form
            ! evaluation based on a generic callback function
            allocate (rtask%p_rfparser)
            call fparser_create(rtask%p_rfparser, itermCount)
            
            ! Fill function parser with data
            istart = 1
            do j = 1,itermCount
              call sys_getNextToken(sparameter, stoken, istart, ',', .false.)
              call fparser_parseFunction(rtask%p_rfparser, j, trim(stoken), (/'x','y','z'/))
            end do
            
            ! That`s it
            exit
          end if
        end do
        
      else
        ! Use unity as constand coefficient
        rtask%p_rform%Dcoefficients(1:rtask%p_rform%itermCount) = 1.0_DP
      end if

      ! Assemble vector data from linear form
      rcollection%p_rfparserQuickAccess1 => rtask%p_rfparser
      if (associated(p_rcubInfo)) then
        call linf_buildVectorScalar(rtask%p_rform, .true.,&
            rtask%p_rvectorScalarDest, p_rcubInfo,&
            problem_buildVectorSc_fparser_sim, rcollection)
      elseif (associated(rtask%p_rspatialDiscr)) then
        call linf_buildVectorScalar(rtask%p_rspatialDiscr,&
            rtask%p_rform, .true., rtask%p_rvectorScalarDest,&
            problem_buildVectorSc_fparser_sim, rcollection)
      else
        call output_line('Neither cubature info structure nor discretisation available',&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
        call sys_halt()
      end if
      
    end subroutine assembleScalarVectorLinf

  end subroutine problem_updateVectorScalar

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateVecBl(rproblemLevel,&
      rvector, rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given block vector on the given problem level with the values
    ! provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which scalar vectors
    ! need to be created and copied
    type(t_vectorTask), pointer :: p_rtasklist

    ! block vector
    type(t_vectorBlock), intent(inout), target :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvectorTemp 
    type(t_vectorTask), pointer :: rtask
    type(t_problem), pointer :: p_rproblemTopLevel
    type(t_problemLevel), pointer :: p_rproblemLevel
    character(len=PARLST_MLDATA) :: skeyword,sparameter,sproblem,stoken,svalue
    character(len=SYS_STRLEN) :: ssubsectionName
    character(len=SYS_STRLEN) :: saction,sperform,sperf
    integer :: iistart,ilev,ivector,istart,itoken,ntoken,imethod,nmethod
    integer(i32) :: iperform

    ! Set pointer to problem structure
    p_rproblemTopLevel => rproblemLevel%p_rproblem
    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec

    ! Get type of action
    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
    call sys_toupper(saction)
    call sys_toupper(sperform)

    ! What action should be performed?
    if (trim(saction) .eq. 'NONE') then
      !-------------------------------------------------------------------------
      ! Do nothing
      
    elseif (trim(saction) .eq. 'CREATE') then

      !-------------------------------------------------------------------------
      ! Create block vector
      !
      ! SYNTAX: vectorblock(n) =
      !           VectorBlock1
      !                ...
      !           VectorBlock2

      ! Create new task
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%bisVectorScalar     = .false.
      rtask%ctask               =  PROBACTION_CREATE
      rtask%ssectionName        =  ssectionName
      rtask%p_rvectorBlockDest  => rvector
      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelDest => rproblemLevel
      
      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this block vector
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
          call sys_halt()
        end if
      end do

      ! Get block discretisation for test and trial space
      call parlst_getvalue_string(rparlist, ssectionName,&
          'discretisation', sparameter)
      call sys_toupper(sparameter)
      rtask%p_rblockDiscr => problem_getBlockDiscr(rproblemLevel%p_rproblem,&
          rproblemLevel, sparameter)
      
      ! Should we perform this task now?
      if (iand(rtask%iperform, iperform) .ne. 0) then
        ! Do we have block discretisations?
        if (.not.associated(rtask%p_rblockDiscr)) then
          call output_line('Unable to find block discretisation',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
          call sys_halt()
        end if
        
        ! Clear preexisting block vector
        call lsysbl_releaseVector(rvector)
        
        ! Create block vector from discretisation
        call lsysbl_createVecBlockByDiscr(rtask%p_rblockDiscr, rvector)
      end if

      !-------------------------------------------------------------------------
      ! Update subvectors
      !
      ! SYNTAX: subvector(n) =
      !           vectorscalar@(i): VectorSc1
      !           vectorscalar@(j): VectorSc2
      !           ...
      !           vectorblock@(a):  VectorBl1
      !
      ! For scalar vectors, the vector objects VectorSc1 and
      ! VectorSc2 are generated at position (i) and (j) of the
      ! given block vector, respectively. For block vectors, the
      ! (1) block of object VectorBl1 is generated at position (a)
      ! of the given vector.  The user is responsible to make sure
      ! that the given vector provides all positions for the subvectors.

      call updateSubvector(rproblemLevel, rvector, rparlist,&
          ssectionName, 0, iperform, p_rtaskList, rproblemTopLevel)
            
    elseif (trim(saction) .eq. 'DUPLICATE') then

      !-------------------------------------------------------------------------
      ! Duplicate block vector are create 1-block wrapper for a scalar vector
      !
      ! SYNTAX: vectorscalar = name,ivector:#,ilev:#,...
      !     OR  vectorblock  = name,ivector:#,ilev:#,...
      !
      ! If ilev is not given then the level of the current problem
      ! level structure is adopted.

      ! Find problem name from which to duplicate scalar vector (if any)
      call parlst_getvalue_string(rparlist, ssectionName, 'vectorscalar', sparameter, '')
      call sys_toupper(sparameter)
      
      if (trim(sparameter) .ne. '') then
        
        !--- Scalar vector case ------------------------------------------------

        ! Create new task for this scalar vector
        nullify(rtask)
        allocate(rtask)
        nullify(rtask%p_rnextTask)
        
        ! Initialise optional parameters
        sproblem = trim(rproblemLevel%p_rproblem%cproblem)
        ilev     = rproblemLevel%ilev
        ivector  = 1
        
        ! Get optional parameters if available
        call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
        
        do itoken = 1,ntoken
          call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
          
          iistart = 1
          call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
          call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
          
          if (trim(skeyword) .eq. 'PROBLEM') then
            sproblem = trim(svalue)
          elseif (trim(skeyword) .eq. 'ILEV') then
            read(svalue,'(I10)') ilev
          elseif (trim(skeyword) .eq. 'IVECTOR') then
            read(svalue,'(I10)') ivector
          else
            call output_line('Unsupported keyword: '//trim(skeyword),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
            call sys_halt()
          end if
        end do
        
        ! Find problem level in global problem structure
        p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
        if (.not.associated(p_rproblemLevel)) then
          call output_line('Unable to find problem level in problem '//trim(sproblem),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
          call sys_halt()
        end if
        
        ! Get section name of scalar vector
        call parlst_getvalue_string(rparlist, trim(sproblem),&
            'vectorscalar', ssubsectionName, isubstring=ivector)
        
        ! Proceed with possibly prerequisite tasks
        call problem_updateVectorScalar(p_rproblemLevel,&
            p_rproblemLevel%RvectorScalar(ivector), rparlist,&
            ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
        
        ! Specify new task for this block vector
        rtask%bisVectorScalar     = .false.
        rtask%ctask               =  PROBACTION_DUPLICATE
        rtask%ssectionName        =  ssectionName
        rtask%p_rvectorScalarSrc  => p_rproblemLevel%RvectorScalar(ivector)
        rtask%p_rvectorBlockDest  => rvector
        rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
        rtask%p_rproblemDest      => rproblemLevel%p_rproblem
        rtask%p_rproblemLevelSrc  => p_rproblemLevel
        rtask%p_rproblemLevelDest => rproblemLevel

      else

        !--- Block vector case -------------------------------------------------

        ! Find problem name from which to duplicate entry of block vector (if any)
        call parlst_getvalue_string(rparlist, ssectionName, 'vectorblock', sparameter, '')
        call sys_toupper(sparameter)
        
        if (trim(sparameter) .ne. '') then
          
          ! Create new task for this scalar vector
          nullify(rtask)
          allocate(rtask)
          nullify(rtask%p_rnextTask)
          
          ! Initialise optional parameters
          sproblem = trim(rproblemLevel%p_rproblem%cproblem)
          ilev     = rproblemLevel%ilev
          ivector  = 1

          ! Get optional parameters if available
          call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
          
          do itoken = 1,ntoken
            call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
            
            iistart = 1
            call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
            call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
            
            if (trim(skeyword) .eq. 'PROBLEM') then
              sproblem = trim(svalue)
            elseif (trim(skeyword) .eq. 'ILEV') then
              read(svalue,'(I10)') ilev
            elseif (trim(skeyword) .eq. 'IVECTOR') then
              read(svalue,'(I10)') ivector
            else
              call output_line('Unsupported keyword: '//trim(skeyword),&
                  OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
              call sys_halt()
            end if
          end do

          ! Find problem level in global problem structure
          p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
          if (.not.associated(p_rproblemLevel)) then
            call output_line('Unable to find problem level in problem '//trim(sproblem),&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
            call sys_halt()
          end if
          
          ! Get section name of block vector
          call parlst_getvalue_string(rparlist, trim(sproblem),&
              'vectorblock', ssubsectionName, isubstring=ivector)
          
          ! Proceed with possibly prerequisite tasks
          call problem_updateVecBl(p_rproblemLevel,&
              p_rproblemLevel%RvectorBlock(ivector), rparlist,&
              ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)

          ! Specify new task for this block vector
          rtask%bisVectorScalar     = .false.
          rtask%ctask               =  PROBACTION_DUPLICATE
          rtask%ssectionName        =  ssectionName
          rtask%p_rvectorBlockSrc   => p_rproblemLevel%RvectorBlock(ivector)
          rtask%p_rvectorBlockDest  => rvector
          rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
          rtask%p_rproblemDest      => rproblemLevel%p_rproblem
          rtask%p_rproblemLevelSrc  => p_rproblemLevel
          rtask%p_rproblemLevelDest => rproblemLevel

        else
          call output_line('Either vectorscalar of vectorblock must be present',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
          call sys_halt()
        end if
      end if
    
      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)

          ! That`s it we do not have to create this block vector
          return
        end if
      else
        p_rtasklist => rtask
      end if      

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
          call sys_halt()
        end if
      end do

      ! Get optional parameters
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupstructure', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupStructure = get_dupType(sparameter)
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sdupcontent', sparameter, 'SHARE')
      call sys_toupper(sparameter)
      rtask%cdupContent = get_dupType(sparameter)

      ! Should we perform this task now?
      if (iand(rtask%iperform, iperform) .ne. 0) then
        if (associated(rtask%p_rvectorScalarSrc)) then
          ! Duplicate scalar vector as 1-block vector
          call lsysbl_createVecFromScalar(rtask%p_rvectorScalarSrc, rvectorTemp)
          
          ! Copy content of temporal vector by hand
          rvector%NEQ         = rvectorTemp%NEQ
          rvector%nblocks     = 1

          allocate(rvector%RvectorBlock(1))
          rvector%RvectorBlock(1) = rvectorTemp%RvectorBlock(1)

          ! Release temporal vector
          call lsysbl_releaseVector(rvectorTemp)

        elseif(associated(rtask%p_rvectorBlockSrc)) then
          ! Duplicate block vector
          call lsysbl_duplicateVector(rtask%p_rvectorBlockSrc,&
              rvector, rtask%cdupStructure, rtask%cdupContent)
        else
          call output_line('Neither scalar nor block vector can be duplicated',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
          call sys_halt()
        end if
      end if
       
    else
      call output_line('Unsupported action: '//trim(saction),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
      call sys_halt()
    end if

    !---------------------------------------------------------------------------
    ! Second part: generate content
    !---------------------------------------------------------------------------
    
    nmethod = parlst_querysubstrings(rparlist, ssectionName, 'smethod')
    allocate(rtask%Cmethod(0:nmethod))

    ! What creation methods are adopted?
    do imethod = 0,nmethod
      
      call parlst_getvalue_string(rparlist, ssectionName,&
          'smethod', sparameter, '', isubstring=imethod)
      call sys_toupper(sparameter)

      if (trim(sparameter) .eq. 'VIRTUAL') then
        
        !-----------------------------------------------------------------------
        ! Create virtual vector, aka, do nothing
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_VIRTUAL

      elseif (trim(sparameter) .eq. 'EMPTY') then
        
        !-----------------------------------------------------------------------
        ! Create empty vector
        
        call lsysbl_clearVector(rvector, 0.0_DP)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_EMPTY
        
      elseif (trim(sparameter) .eq. 'UNITY') then
        
        !-----------------------------------------------------------------------
        ! Create unit vector
        
        call lsysbl_clearVector(rvector, 1.0_DP)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_UNITY
        
      elseif (trim(sparameter) .eq. 'DOF') then
        
        !-----------------------------------------------------------------------
        ! Create coordinate vector for degrees of freedom
        
        call lin_calcDofCoordsBlock(rtask%p_rblockDiscr, rvector)
        
        ! Update internal data of the task item
        rtask%Cmethod(imethod) = PROBMETH_VECTOR_DOF
        
      else
        call output_line('Unsupported method: '//trim(sparameter),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
        call sys_halt()
      end if
    end do

  contains

    ! Here, some auxiliary routines follow

    !***************************************************************************
    ! Search for given task in the list of tasks. If the task is not
    ! present in the list, then it is appended.
    function appendToTaskList(rtasklist, rtask) result(bexists)

      type(t_vectorTask), intent(inout), target :: rtasklist
      type(t_vectorTask), intent(in), target :: rtask
      
      ! local variable
      type(t_vectorTask), pointer :: p_rtask
      logical :: bexists      

      p_rtask => rtasklist
      do while(associated(p_rtask))
        if (associated(p_rtask%p_rvectorBlockDest,&
            rtask%p_rvectorBlockDest)) then
          bexists = .true.
          return
        end if
        if (.not.associated(p_rtask%p_rnextTask)) then
          p_rtask%p_rnextTask => rtask
          bexists = .false.
          return
        else
          p_rtask => p_rtask%p_rnextTask
        end if
      end do
      
    end function appendToTaskList

    !**************************************************************
    ! Return the numeric value of the duplication flad
    function get_dupType(sdupType) result(cdupType)

      character(len=*), intent(in) :: sdupType
      integer :: cdupType

      if (sdupType .eq. 'IGNORE') then
        cdupType = LSYSSC_DUP_IGNORE
      elseif (sdupType .eq. 'REMOVE') then
        cdupType = LSYSSC_DUP_REMOVE
      elseif (sdupType .eq. 'DISMISS') then
        cdupType = LSYSSC_DUP_DISMISS
      elseif (sdupType .eq. 'SHARE') then
        cdupType = LSYSSC_DUP_SHARE
      elseif (sdupType .eq. 'COPY') then
        cdupType = LSYSSC_DUP_COPY
      elseif (sdupType .eq. 'COPYOVERWRITE') then
        cdupType = LSYSSC_DUP_COPYOVERWRITE
      elseif (sdupType .eq. 'ASIS') then
        cdupType = LSYSSC_DUP_ASIS
      elseif (sdupType .eq. 'EMPTY') then
        cdupType = LSYSSC_DUP_EMPTY
      elseif (sdupType .eq. 'TEMPLATE') then
        cdupType = LSYSSC_DUP_TEMPLATE
      else
        cdupType = 0
      end if

    end function get_dupType
    
    !**************************************************************
    ! Update the subvector
    recursive subroutine updateSubvector(rproblemLevel, rvector,&
        rparlist, ssectionName, iblockOffset, iperform,&
        p_rtaskList, rproblemTopLevel)

      type(t_problemLevel), intent(in), target :: rproblemLevel
      type(t_problem), intent(in), target, optional :: rproblemTopLevel
      type(t_vectorBlock), intent(inout) :: rvector
      type(t_parlist), intent(in) :: rparlist
      character(len=*), intent(in) :: ssectionName
      integer, intent(in) :: iblockOffset
      integer(i32), intent(in) :: iperform
      type(t_vectorTask), pointer :: p_rtaskList

      ! local variables
      character(len=PARLST_MLDATA) :: sparameter,stoken,ssubtoken,ssubvector
      integer :: i,iblock,nsubvector,istart,iistart,iiistart

      ! Get number of subvectors
      nsubvector = parlst_querysubstrings(rparlist, ssectionName, 'subvector')

      ! Loop over all subvectors
      do i=1,nsubvector
        
        call parlst_getvalue_string(rparlist, ssectionName, 'subvector',&
            sparameter, isubstring=i)
        call sys_toupper(sparameter)

        ! Extract content from parameter file line by line
        ! SYNTAX: vectorscalar@(1): VectorSc1
        !     OR  vectorblock@(1):  VectorBl1
        istart = 1; iistart = 1; iiistart = 1
        call sys_getNextToken (sparameter, stoken, istart, ":", .false.)
        call sys_getNextToken (stoken, ssubvector, iistart, "@", .false.)
        call sys_getNextToken (stoken, ssubtoken, iiistart, ",", .false.)
        read(ssubtoken(iistart:iiistart),*) iblock
        iistart = iiistart
        call sys_getNextToken (stoken, ssubtoken, iiistart, ")", .false.)
        read(ssubtoken,*) iblock
        call sys_getNextToken (sparameter, stoken, istart, ":", .false.)

        ! Apply given block offsets
        iblock = iblock+iblockOffset

        ! Check of this is a scalar vector or a block vector
        if (trim(ssubvector) .eq. "VECTORSCALAR") then
          
          ! Check if scalar subvector is available
          if ((rvector%nblocks .lt. iblock)) then
            call output_line('Position of subvector is not available',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
            call sys_halt()
          end if

          ! Update single scalar subvector
          call problem_updateVectorScalar(rproblemLevel,&
              rvector%RvectorBlock(iblock), rparlist, trim(stoken),&
              p_rtasklist, rproblemTopLevel, iperform)
          
        elseif(trim(ssubvector) .eq. "VECTORBLOCK") then
          
          ! Recursively update block vector
          call updateSubvector(rproblemLevel, rvector, rparlist,&
              trim(stoken), iblock-1, iperform, p_rtasklist, rproblemTopLevel)
        else
          call output_line('Either vectorscalar of vectorblock must be present',&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVecBl')
          call sys_halt()
        end if
        
      end do

    end subroutine updateSubvector

  end subroutine problem_updateVecBl

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_updateVector(rtasklist, iperformSpec)

!<description>
    ! This subroutine updates the scalar/block vectors with the
    ! internal values assigned to the items of the task list.
!</description>

!<input>
    ! task list
    type(t_vectorTask), intent(in), target :: rtasklist

    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>
!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type(t_vectorBlock) :: rvectorTemp
    type(t_vectorTask), pointer :: p_rtask
    integer :: iperform,imethod

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec
    
    ! Iterate over all tasks
    p_rtask => rtasklist
    taskloop: do while (associated(p_rtask))

      ! Do we have to perform this task?
      if (iand(p_rtask%iperform, iperform) .eq. 0) then
        p_rtask => p_rtask%p_rnextTask
        cycle taskloop
      end if
      
      !-------------------------------------------------------------------------
      ! First part: create or duplicate vector
      !-------------------------------------------------------------------------
      select case(p_rtask%ctask)
      case(PROBACTION_CREATE)

        ! Do we have to assemble a scalar vector?
        if (p_rtask%bisVectorScalar) then
          
          !---------------------------------------------------------------------
          ! Create scalar vector
          
          ! Do we have spatial discretisations?
          if (.not.associated(p_rtask%p_rspatialDiscr)) then
            call output_line('Unable to find spatial discretisation',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
            call sys_halt()
          end if
          
          ! Clear preexisting vector
          call lsyssc_releaseVector(p_rtask%p_rvectorScalarDest)
          
          ! Create vector structure and set data type
          call lsyssc_createVector(p_rtask%p_rspatialDiscr,&
              p_rtask%p_rvectorScalarDest, p_rtask%nvar, .false., p_rtask%cdataType)

        else

          !---------------------------------------------------------------------
          ! Create block vector

          ! Do we have block discretisations?
          if (.not.associated(p_rtask%p_rblockDiscr)) then
            call output_line('Unable to find block discretisation',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
            call sys_halt()
          end if

          ! Clear preexisting block vector
          call lsysbl_releaseVector(p_rtask%p_rvectorBlockDest)

          ! Create block vector from discretisation
          call lsysbl_createVecBlockByDiscr(p_rtask%p_rblockDiscr,&
              p_rtask%p_rvectorBlockDest)
          
          ! Proceed with next task
          p_rtask => p_rtask%p_rnextTask
          
          ! That`s it for block vectors
          cycle taskloop
        end if

      case(PROBACTION_DUPLICATE)

        !-----------------------------------------------------------------------
        ! Duplicate scalar/block vector

        if (p_rtask%bisVectorScalar) then
          ! Duplicate scalar vector
          call lsyssc_duplicateVector(p_rtask%p_rvectorScalarSrc,&
              p_rtask%p_rvectorScalarDest, p_rtask%cdupStructure, p_rtask%cdupContent)
        else
          if (associated(p_rtask%p_rvectorScalarSrc)) then
            ! Duplicate scalar vector as 1-block vector
            call lsysbl_createVecFromScalar(p_rtask%p_rvectorScalarSrc, rvectorTemp)

            ! Copy content of temporal vector by hand
            p_rtask%p_rvectorBlockDest%NEQ         = rvectorTemp%NEQ
            p_rtask%p_rvectorBlockDest%nblocks     = 1

            allocate(p_rtask%p_rvectorBlockDest%RvectorBlock(1))
            p_rtask%p_rvectorBlockDest%RvectorBlock(1) = rvectorTemp%RvectorBlock(1)

            ! Release temporal vector
            call lsysbl_releaseVector(rvectorTemp)

          elseif(associated(p_rtask%p_rvectorBlockSrc)) then
            ! Duplicate block vector
            call lsysbl_duplicateVector(p_rtask%p_rvectorBlockSrc,&
                p_rtask%p_rvectorBlockDest, p_rtask%cdupStructure, p_rtask%cdupContent)
          else
            call output_line('Neither scalar nor block vector can be duplicated',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
            call sys_halt()
          end if
        end if

      case default
        call output_line('Unsupported action: '//sys_siL(p_rtask%ctask,3),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMamtrix')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------
      ! Second part: generate content (only for scalar vectors)
      !-------------------------------------------------------------------------

       ! What creation methods are adopted?
      do imethod = lbound(p_rtask%Cmethod,1), ubound(p_rtask%Cmethod,1)
        
        select case (p_rtask%Cmethod(imethod))
        case(PROBMETH_VECTOR_EMPTY)
          
          !-------------------------------------------------------------------
          ! Create empty vector
          
          call lsyssc_clearVector(p_rtask%p_rvectorScalarDest, 0.0_DP)
          
        case(PROBMETH_VECTOR_UNITY)
          
          !-------------------------------------------------------------------
          ! Create unit vector
          
          call lsyssc_clearVector(p_rtask%p_rvectorScalarDest, 1.0_DP)
          
        case(PROBMETH_VECTOR_LINF)
          
          !-------------------------------------------------------------------
          ! Create vector by bilinearform evaluation         
          
          rcollection%p_rfparserQuickAccess1 => p_rtask%p_rfparser
          if (associated(p_rtask%p_rcubatureInfo)) then
            call linf_buildVectorScalar(p_rtask%p_rform, .true.,&
                p_rtask%p_rvectorScalarDest, p_rtask%p_rcubatureInfo,&
                problem_buildVectorSc_fparser_sim, rcollection)
          elseif (associated(p_rtask%p_rspatialDiscr)) then
            call linf_buildVectorScalar(p_rtask%p_rspatialDiscr,&
                p_rtask%p_rform, .true., p_rtask%p_rvectorScalarDest,&
                problem_buildVectorSc_fparser_sim, rcollection)
          else
            call output_line('Neither cubature info structure nor discretisation available',&
                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
            call sys_halt()
          end if
          
        case(PROBMETH_VECTOR_VIRTUAL)
          
          !-------------------------------------------------------------------
          ! Create virtual vector, aka, do nothing
          
        case default
          call output_line('Unsupported method: '//sys_siL(p_rtask%Cmethod(imethod),3),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
          call sys_halt()
        end select
      end do

      ! Proceed with next task
      p_rtask => p_rtask%p_rnextTask
    end do taskloop
    
  end subroutine problem_updateVector

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_clearVectorTL(p_rtasklist)

!<description>
    ! This subroutine clears the scalar/block vector task list.
!</description>

!<input>
    ! task list
    type(t_vectorTask), pointer :: p_rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_vectorTask), pointer :: p_rtask,p_rtaskSave

    p_rtask => p_rtasklist
    do while(associated(p_rtask))
      if (associated(p_rtask%Cmethod)) deallocate(p_rtask%Cmethod)
      if (associated(p_rtask%p_rform)) deallocate(p_rtask%p_rform)
      if (associated(p_rtask%p_rfparser)) then
        call fparser_release(p_rtask%p_rfparser)
        deallocate(p_rtask%p_rfparser)
      end if
      p_rtaskSave => p_rtask
      p_rtask     => p_rtask%p_rnextTask
      deallocate(p_rtaskSave)
    end do
    
  end subroutine problem_clearVectorTL

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_infoVectorTL(rtasklist)

!<description>
    ! This subroutine outputs information about the scalar/block
    ! vector task list.
!</description>

!<input>
    ! task list
    type(t_vectorTask), intent(in), target :: rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_vectorTask), pointer :: p_rtask
    integer :: i

    p_rtask => rtasklist
    do while(associated(p_rtask))

      select case(p_rtask%ctask)
      case (PROBACTION_NONE)
        call output_lbrk()
        call output_line ('VectorTask: DO NOTHING')

      case (PROBACTION_CREATE)
        call output_lbrk()
        call output_line ('VectorTask: CREATE')

      case (PROBACTION_DUPLICATE)
        call output_lbrk()
        call output_line ('VectorTask: DUPLICATE')

      case default
        call output_line ('VectorTask: UNSUPPORTED')
      end select

      call output_line ('-----------------------')
      call output_line ('vector type               : '//merge('SCALAR', 'BLOCK ', p_rtask%bisVectorScalar))
      call output_line ('section name              : '//trim(p_rtask%ssectionName))
      call output_line ('destination problem       : '//trim(p_rtask%p_rproblemDest%cproblem))
      call output_line ('destination problem level : '//trim(sys_siL(p_rtask%p_rproblemLevelDest%ilev,15)))
      call output_line ('destination vector scalar : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rvectorScalarDest)))
      call output_line ('destination vector block  : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rvectorBlockDest)))
      if (associated(p_rtask%p_rproblemSrc))&
      call output_line ('source problem            : '//trim(p_rtask%p_rproblemSrc%cproblem))
      if (associated(p_rtask%p_rproblemLevelSrc))&
      call output_line ('source problem level      : '//trim(sys_siL(p_rtask%p_rproblemLevelSrc%ilev,15)))
      call output_line ('source vector scalar      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rvectorScalarSrc)))
      call output_line ('source vector block       : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rvectorBlockSrc)))
      call output_line ('iperform                  : '//trim(sys_siL(int(p_rtask%iperform),15)))
      call output_line ('method                    : '//trim(sys_siL(p_rtask%cmethod(0),15)))
      do i = 1,ubound(p_rtask%Cmethod,1)
        call output_line ('                          : '//trim(sys_siL(p_rtask%Cmethod(i),15)))
      end do
      call output_line ('data type                 : '//trim(sys_siL(p_rtask%cdatatype,15)))
      call output_line ('nvar                      : '//trim(sys_siL(p_rtask%nvar,15)))
      call output_line ('spatial discretisation:     '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rspatialDiscr)))
      call output_line ('block discretisation      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rblockDiscr)))
      call output_line ('cubature info structure   : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rcubatureInfo)))
      call output_line ('function parser           : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rfparser)))
      call output_line ('linear form evaluator     : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
                                                              associated(p_rtask%p_rform)))

      p_rtask => p_rtask%p_rnextTask
    end do
    
  end subroutine problem_infoVectorTL

  !*****************************************************************************
  ! AUXILIARY SUBROUTINES
  !*****************************************************************************

!<function>
  function problem_getSpatialDiscr(rproblem, rproblemLevel, sspatialDiscr,& 
      rspatialDiscrDefault) result(p_rspatialDiscr)

!<descrition>
    ! This subroutines sets the pointer to a spatial discretisation
    ! structure based on the string representation sspatialDiscr and the
    ! data present in the problem structure
!</descrition>
   
!<input>
    ! problem structure
    type(t_problem), intent(in) :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! string representation
    character(len=*), intent(in) :: sspatialDiscr

    ! OPTIONAL: default spatial discretisation structure to be used if
    ! the specified structure cannot be found
    type(t_spatialDiscretisation), intent(in), target, optional :: rspatialDiscrDefault
!</input>

!<result>
    ! pointer to the specified spatial discretisation structure
    ! or null of it is not available
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
!</result>

!</function>

      ! local variables
      type(t_problemLevel), pointer :: p_rproblemLevel
      type(t_blockDiscretisation), pointer :: p_rblockDiscr
      character(len=PARLST_MLDATA) :: skeyword,sproblem,stoken,svalue
      integer :: iblock,idiscr,iistart,ilev,istart,itoken,ntoken

      ! Do we have a valid string?
      if (trim(sspatialDiscr) .eq. '') then
        if (present(rspatialDiscrDefault)) then
          p_rspatialDiscr => rspatialDiscrDefault
        else
          nullify(p_rspatialDiscr)
        end if
        ! That`s it
      end if

      ! Initialise parameters
      sproblem = trim(rproblemLevel%p_rproblem%cproblem)
      ilev     = rproblemLevel%ilev
      idiscr   = 1
      iblock   = 1
      
      call sys_countTokens(sspatialDiscr, ntoken, ',', .false.); istart=1
      
      do itoken = 1,ntoken
        call sys_getNextToken (sspatialDiscr, stoken, istart, ',', .false.)
        
        iistart = 1
        call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
        call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
        
        if (trim(skeyword) .eq. 'PROBLEM') then
          sproblem = trim(svalue)
        elseif (trim(skeyword) .eq. 'ILEV') then
          read(svalue,'(I10)') ilev
        elseif (trim(skeyword) .eq. 'IBLOCK') then
          read(svalue,'(I10)') iblock
        elseif (trim(skeyword) .eq. 'IDISCR' .or.&
            trim(skeyword) .eq. 'IDISCRETISATION') then
          read(svalue,'(I10)') idiscr
        else
          call output_line('Unsupported keyword: '//trim(skeyword),&
              OU_CLASS_ERROR,OU_MODE_STD,'poblem_getSpatialDiscr')
          call sys_halt()
        end if
      end do
      
      ! Find problem level in global problem structure
      p_rproblemLevel => problem_getLevel(rproblem, trim(sproblem), ilev)       
      if (.not.associated(p_rproblemLevel)) then
        call output_line('Unable to find problem level in problem '//trim(sproblem),&
            OU_CLASS_ERROR,OU_MODE_STD,'poblem_getSpatialDiscr')
        call sys_halt()
      end if
      
      ! Set pointer to discretisation structure
      nullify(p_rspatialDiscr)
      if ((idiscr .ge. lbound(p_rproblemLevel%Rdiscretisation,1)) .and.&
          (idiscr .le. ubound(p_rproblemLevel%Rdiscretisation,1))) then
        p_rblockDiscr => p_rproblemLevel%Rdiscretisation(idiscr)
        if ((iblock .ge. lbound(p_rblockDiscr%RspatialDiscr,1)) .and.&
            (iblock .le. ubound(p_rblockDiscr%RspatialDiscr,1))) then
          p_rspatialDiscr => p_rblockDiscr%RspatialDiscr(iblock)
        end if
      end if
      
    end function problem_getSpatialDiscr

    !*****************************************************************************

!<function>
  function problem_getBlockDiscr(rproblem, rproblemLevel, sblockDiscr,& 
      rblockDiscrDefault) result(p_rblockDiscr)

!<descrition>
    ! This subroutines sets the pointer to a block discretisation
    ! structure based on the string representation sblockDiscr and the
    ! data present in the problem structure
!</descrition>
   
!<input>
    ! problem structure
    type(t_problem), intent(in) :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! string representation
    character(len=*), intent(in) :: sblockDiscr

    ! OPTIONAL: default block discretisation structure to be used if
    ! the specified structure cannot be found
    type(t_blockDiscretisation), intent(in), target, optional :: rblockDiscrDefault
!</input>

!<result>
    ! pointer to the specified block discretisation structure
    ! or null of it is not available
    type(t_blockDiscretisation), pointer :: p_rblockDiscr
!</result>

!</function>

      ! local variables
      type(t_problemLevel), pointer :: p_rproblemLevel
      character(len=PARLST_MLDATA) :: skeyword,sproblem,stoken,svalue
      integer :: idiscr,iistart,ilev,istart,itoken,ntoken

      ! Do we have a valid string?
      if (trim(sblockDiscr) .eq. '') then
        if (present(rblockDiscrDefault)) then
          p_rblockDiscr => rblockDiscrDefault
        else
          nullify(p_rblockDiscr)
        end if
        ! That`s it
      end if

      ! Initialise parameters
      sproblem = trim(rproblemLevel%p_rproblem%cproblem)
      ilev     = rproblemLevel%ilev
      idiscr   = 1
      
      call sys_countTokens(sblockDiscr, ntoken, ',', .false.); istart=1
      
      do itoken = 1,ntoken
        call sys_getNextToken (sblockDiscr, stoken, istart, ',', .false.)
        
        iistart = 1
        call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
        call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
        
        if (trim(skeyword) .eq. 'PROBLEM') then
          sproblem = trim(svalue)
        elseif (trim(skeyword) .eq. 'ILEV') then
          read(svalue,'(I10)') ilev
        elseif (trim(skeyword) .eq. 'IDISCR' .or.&
            trim(skeyword) .eq. 'IDISCRETISATION') then
          read(svalue,'(I10)') idiscr
        else
          call output_line('Unsupported keyword: '//trim(skeyword),&
              OU_CLASS_ERROR,OU_MODE_STD,'poblem_getBlockDiscr')
          call sys_halt()
        end if
      end do
      
      ! Find problem level in global problem structure
      p_rproblemLevel => problem_getLevel(rproblem, trim(sproblem), ilev)       
      if (.not.associated(p_rproblemLevel)) then
        call output_line('Unable to find problem level in problem '//trim(sproblem),&
            OU_CLASS_ERROR,OU_MODE_STD,'poblem_getBlockDiscr')
        call sys_halt()
      end if
      
      ! Set pointer to discretisation structure
      nullify(p_rblockDiscr)
      if ((idiscr .ge. lbound(p_rproblemLevel%Rdiscretisation,1)) .and.&
          (idiscr .le. ubound(p_rproblemLevel%Rdiscretisation,1))) then
        p_rblockDiscr => p_rproblemLevel%Rdiscretisation(idiscr)        
      end if
      
    end function problem_getBlockDiscr

    !*****************************************************************************

!<function>
  function problem_getCubInfo(rproblem, rproblemLevel, scubInfo,& 
      rcubInfoDefault) result(p_rcubInfo)

!<descrition>
    ! This subroutines sets the pointer to a cubature info structure
    !  based on the string representation scubInfo and the data
    !  present in the problem structure
!</descrition>
   
!<input>
    ! problem structure
    type(t_problem), intent(in) :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! string representation
    character(len=*), intent(in) :: scubInfo

    ! OPTIONAL: default cubature info structure to be used if
    ! the specified structure cannot be found
    type(t_scalarCubatureInfo), intent(in), target, optional :: rcubInfoDefault
!</input>

!<result>
    ! pointer to the specified cubature info structure
    ! or null of it is not available
    type(t_scalarCubatureInfo), pointer :: p_rcubInfo
!</result>

!</function>

    ! local variables
    type(t_problemLevel), pointer :: p_rproblemLevel
    character(len=PARLST_MLDATA) :: skeyword,sproblem,stoken,svalue
    integer :: icubinfo,iistart,ilev,istart,itoken,ntoken
    
    ! Do we have a valid string?
    if (trim(scubInfo) .eq. '') then
      if (present(rcubInfoDefault)) then
        p_rcubInfo => rcubInfoDefault
      else
        nullify(p_rcubInfo)
      end if
      ! That`s it
      return
    end if
    
    ! Initialise parameters
    sproblem = trim(rproblemLevel%p_rproblem%cproblem)
    ilev     = rproblemLevel%ilev
    icubinfo = 1
    
    call sys_countTokens(scubInfo, ntoken, ',', .false.); istart=1
    
    do itoken = 1,ntoken
      call sys_getNextToken (scubInfo, stoken, istart, ',', .false.)
      
      iistart = 1
      call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
      call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
      
      if (trim(skeyword) .eq. 'PROBLEM') then
        sproblem = trim(svalue)
      elseif (trim(skeyword) .eq. 'ILEV') then   
        read(svalue,'(I10)') ilev
      elseif (trim(skeyword) .eq. 'ICUBINFO' .or.&
          trim(skeyword) .eq. 'ICUBATUREINFO') then
        read(svalue,'(I10)') icubinfo
      else
        call output_line('Unsupported keyword: '//trim(skeyword),&
            OU_CLASS_ERROR,OU_MODE_STD,'problem_getCubInfo')
        call sys_halt()
      end if
    end do
    
    ! Find problem level in global problem structure
    p_rproblemLevel => problem_getLevel(rproblem, trim(sproblem), ilev)       
    if (.not.associated(p_rproblemLevel)) then
      call output_line('Unable to find problem level in problem '//trim(sproblem),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_getCubInfo')
      call sys_halt()
    end if

    ! Do we have a cubature info?
    if ((icubinfo .ge. lbound(p_rproblemLevel%Rcubatureinfo,1)) .and.&
        (icubinfo .le. ubound(p_rproblemLevel%Rcubatureinfo,1))) then
      p_rcubInfo => p_rproblemLevel%Rcubatureinfo(icubinfo)
    else
      call output_line('Unable to find cubature info structure',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_getCubInfo')
      call sys_halt()
    end if
    
  end function problem_getCubInfo
  
  !****************************************************************************

!<subroutine>

  subroutine problem_buildMatrixSc_fparser_sim (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, IdofsTrial, IdofsTest, rdomainIntSubset,&
                  Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use spatialdiscretisation
    use triangulation
    use fsystem

  !<description>
    ! This subroutine can be called during the matrix assembly.
    ! It computes the coefficients in front of the terms of the
    ! bilinear form from the given function parser.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>

  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear vectors).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  !</input>

  !<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>

  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>

  !</subroutine>

    ! local variable
    type(t_fparser), pointer :: p_rfparser
    integer :: i,iel,ipoint

    ! Get function parser from collection
    p_rfparser => rcollection%p_rfparserQuickAccess1

    if (associated(p_rfparser)) then
      
      do iel=1,size(Dcoefficients,3)
        do ipoint=1,size(Dcoefficients,2)
          do i=1,size(Dcoefficients,1)
            call fparser_evalFunction(p_rfparser, i, Dpoints(:,ipoint,iel),&
                                      Dcoefficients(i,ipoint,iel))
          end do
        end do
      end do

    else
      call output_line('Quick access function parse is not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_buildMatrixSc_fparser_sim')
      call sys_halt()
    end if

  end subroutine problem_buildMatrixSc_fparser_sim

  !****************************************************************************

!<subroutine>

  subroutine problem_buildVectorSc_fparser_sim (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, &
                  IdofsTest, rdomainIntSubset,&
                  Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use spatialdiscretisation
    use triangulation
    use fsystem

  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>

  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements test in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  !</input>

  !<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>

  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>

  !</subroutine>

    ! local variable
    type(t_fparser), pointer :: p_rfparser
    integer :: i,iel,ipoint

    ! Get function parser from collection
    p_rfparser => rcollection%p_rfparserQuickAccess1

    if (associated(p_rfparser)) then
      
      do iel=1,size(Dcoefficients,3)
        do ipoint=1,size(Dcoefficients,2)
          do i=1,size(Dcoefficients,1)
            call fparser_evalFunction(p_rfparser, i, Dpoints(:,ipoint,iel),&
                                      Dcoefficients(i,ipoint,iel))
          end do
        end do
      end do

    else

      do iel=1,size(Dcoefficients,3)
        do ipoint=1,size(Dcoefficients,2)
          do i=1,size(Dcoefficients,1)
            Dcoefficients(i,ipoint,iel) = rform%Dcoefficients(i)
          end do
        end do
      end do
    end if

  end subroutine problem_buildVectorSc_fparser_sim

  !*****************************************************************************

!<subroutine>

  subroutine problem_initGroupFEMBlockAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises all blocks of group finite element
    ! sets of the given problem level with the values provided by the
    ! parameter list.  If the optional parameter p_rtasklist is given,
    ! then the task list which is created during the initialisation is
    ! returned.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! OPTIONAL: task list
    type(t_groupFEMBlockTask), intent(inout), pointer, optional :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_groupFEMBlockTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update all blocks of group finite element structures
    call problem_updateGroupFEMBlockAll(rproblemLevel,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearGroupFEMBlockTL(p_rtasklistLocal)
    end if

  end subroutine problem_initGroupFEMBlockAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_updateGroupFEMBlockAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine the task list for all blocks of group finite
    ! element sets of the given problem level with the values provided
    ! by the parameter list.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! pointer to the task list
    type(t_groupFEMBlockTask), pointer :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: ssubsectionName
    integer :: i,ngroupfemBlock

    ! Consistency check
    ngroupfemBlock = parlst_querysubstrings(rparlist, ssectionName,&
                                            'groupfemblock')
    if ((ngroupfemBlock .gt. 0) .and.&
        (ngroupfemBlock .ne. size(rproblemLevel%RgroupFEMBlock))) then
      call output_line('Invalid number of blocks of group finite element sets',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateGroupFEMBlockAll')
      call sys_halt()
    end if
   
    ! Update all blocks of group finite element sets
    do i=1,ngroupfemblock
      call parlst_getvalue_string(rparlist, ssectionName,&
          'groupfemblock', ssubsectionName, isubstring=i)
      call problem_updateGroupFEMBlock1(rproblemLevel,&
          rproblemLevel%RgroupFEMBlock(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do
    
  end subroutine problem_updateGroupFEMBlockAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_initGroupFEMBlock(rproblemLevel, rgroupFEMBlock,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given block of group finite
    ! element sets on the given problem level with the values provided
    ! by the parameter list. If the optional parameter p_rtasklist is
    ! given, then the task list which is created during the
    ! initialisation is returned.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which block of group
    ! finite element sets to be created and copied
    type(t_groupFEMBlockTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! block of group finite element sets
    type(t_groupFEMBlock), intent(out), target :: rgroupFEMBlock
!</output>
!</subroutine>

    ! local variables
    type(t_groupFEMBlockTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given block of group finite element sets
    call problem_updateGroupFEMBlock1(rproblemLevel, rgroupFEMBlock,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearGroupFEMBlockTL(p_rtasklistLocal)
    end if

  end subroutine problem_initGroupFEMBlock

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateGroupFEMBlock1(rproblemLevel,&
      rgroupFEMBlock, rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given block of group finite element sets on the given problem
    ! level with the values provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which block of group
    ! finite element sets need to be created and copied
    type(t_groupFEMBlockTask), pointer :: p_rtasklist

    ! @block of group finite element sets
    type(t_groupFEMBlock), intent(inout), target :: rgroupFEMBlock
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvectorBlock
    type(t_groupFEMBlockTask), pointer :: rtask
    type(t_problem), pointer :: p_rproblemTopLevel
    type(t_problemLevel), pointer :: p_rproblemLevel
    character(len=PARLST_MLDATA) :: skeyword,sparameter,sproblem,stoken,svalue
    character(len=SYS_STRLEN) :: ssubsectionName
    character(len=SYS_STRLEN) :: saction,sperform,sperf
    integer :: iblock,iistart,ilev,ivector,istart,itoken,ntoken,nblocks
    integer(i32) :: iperform

    ! Set pointer to problem structure
    p_rproblemTopLevel => rproblemLevel%p_rproblem
    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel

    iperform = PROBACTION_PERFORM_ALWAYS
    if (present(iperformSpec)) iperform = iperformSpec

    ! Get type of action
    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
    call sys_toupper(saction)
    call sys_toupper(sperform)

    !---------------------------------------------------------------------------
    ! First part: generate structure
    !---------------------------------------------------------------------------

    ! What action should be performed?
    if (trim(saction) .eq. 'NONE') then
      !-------------------------------------------------------------------------
      ! Do nothing
      
    elseif (trim(saction) .eq. 'CREATE') then

      !-------------------------------------------------------------------------
      ! Create block of group finite element sets
      !
      ! SYNTAX: groupfemblock(n) =
      !           GroupFEMSet1
      !                ...
      !           GroupFEMSet2

      ! Create new task
      nullify(rtask)
      allocate(rtask)
      nullify(rtask%p_rnextTask)
      rtask%ctask                =  PROBACTION_CREATE
      rtask%ssectionName         =  ssectionName
      rtask%p_rgroupFEMBlockDest => rgroupFEMBlock
      rtask%p_rproblemDest       => rproblemLevel%p_rproblem
      rtask%p_rproblemLevelDest  => rproblemLevel

      ! Append task to task list (if not already present)
      if (associated(p_rtasklist)) then
        if (appendToTaskList(p_rtasklist, rtask)) then
          deallocate(rtask)
          
          ! That`s it we do not have to create this scalar vector
          return
        end if
      else
        p_rtasklist => rtask
      end if

      ! When should we perform this task?
      call sys_countTokens(sperform, ntoken, ',', .false.)
      istart = 1
      do itoken=1,max(ntoken,1)
        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
        if (trim(adjustl(sperf)) .eq. 'INIT') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
        else
          call output_line('Unsupported perform type: '//trim(sperf),&
              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateGroupFEMBlock')
          call sys_halt()
        end if
      end do
!!$      
!!$      ! Get optional parameters
!!$    call parlst_getvalue_int(rparlist, ssectionName,&
!!$        'cdatatype', rtask%cdataType, ST_NOHANDLE)
!!$    if (rtask%cdataType .eq. ST_NOHANDLE) then
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'sdatatype', sparameter, 'DOUBLE')
!!$      call sys_toupper(sparameter)
!!$      rtask%cdataType = get_dataType(trim(sparameter))
!!$    end if
!!$      call parlst_getvalue_int(rparlist, ssectionName,&
!!$          'nvar', rtask%nvar, 1)
!!$      
!!$      ! Get spatial discretisation for test and trial space
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'discretisation', sparameter)
!!$      call sys_toupper(sparameter)
!!$      rtask%p_rspatialDiscr => problem_getSpatialDiscr(&
!!$          rproblemLevel%p_rproblem, rproblemLevel, sparameter)
!!$      
!!$      ! Should we perform this task now?
!!$      if (iand(rtask%iperform, iperform) .ne. 0) then
!!$        ! Do we have spatial discretisations?
!!$        if (.not.associated(rtask%p_rspatialDiscr)) then
!!$          call output_line('Unable to find spatial discretisation',&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$        
!!$        ! Clear preexisting scalar vector
!!$        call lsyssc_releaseVector(rvector)
!!$        
!!$        ! Create vector structure and set data type
!!$        call lsyssc_createVector(rtask%p_rspatialDiscr, rvector,&
!!$            rtask%nvar, .false., rtask%cdataType)
!!$      end if
!!$
!!$    elseif (trim(saction) .eq. 'DUPLICATE') then
!!$
!!$      !-------------------------------------------------------------------------
!!$      ! Duplicate scalar vector or subvector of a block vector
!!$      !
!!$      ! SYNTAX: vectorscalar = name,ivector:#,ilev:#,...
!!$      !     OR  vectorblock  = name,ivector:#,ilev:#,iblock:#
!!$      !
!!$      ! If ilev is not given then the level of the current problem
!!$      ! level structure is adopted. If iblock are not given then
!!$      ! standard value 1 is used in both cases
!!$
!!$      ! Find problem name from which to duplicate scalar vector (if any)
!!$      call parlst_getvalue_string(rparlist, ssectionName, 'vectorscalar', sparameter, '')
!!$      call sys_toupper(sparameter)
!!$
!!$      if (trim(sparameter) .ne. '') then
!!$        
!!$        !--- Scalar vector case ------------------------------------------------
!!$
!!$        ! Create new task for this scalar vector
!!$        nullify(rtask)
!!$        allocate(rtask)
!!$        nullify(rtask%p_rnextTask)
!!$        
!!$        ! Initialise optional parameters
!!$        sproblem = trim(rproblemLevel%p_rproblem%cproblem)
!!$        ilev     = rproblemLevel%ilev
!!$        ivector  = 1
!!$        
!!$        ! Get optional parameters if available
!!$        call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
!!$        
!!$        do itoken = 1,ntoken
!!$          call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
!!$          
!!$          iistart = 1
!!$          call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
!!$          call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
!!$          
!!$          if (trim(skeyword) .eq. 'PROBLEM') then
!!$            sproblem = trim(svalue)
!!$          elseif (trim(skeyword) .eq. 'ILEV') then
!!$            read(svalue,'(I10)') ilev
!!$          elseif (trim(skeyword) .eq. 'IVECTOR') then
!!$            read(svalue,'(I10)') ivector
!!$          else
!!$            call output_line('Unsupported keyword: '//trim(skeyword),&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$            call sys_halt()
!!$          end if
!!$        end do
!!$        
!!$        ! Find problem level in global problem structure
!!$        p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
!!$        if (.not.associated(p_rproblemLevel)) then
!!$          call output_line('Unable to find problem level in problem '//trim(sproblem),&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$        
!!$        ! Get section name of scalar vector
!!$        call parlst_getvalue_string(rparlist, trim(sproblem),&
!!$            'vectorscalar', ssubsectionName, isubstring=ivector)
!!$
!!$        ! Proceed with possibly prerequisite tasks
!!$        call problem_updateVectorScalar(p_rproblemLevel,&
!!$            p_rproblemLevel%RvectorScalar(ivector), rparlist,&
!!$            ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
!!$        
!!$        ! Specify new task for this cubature info structure
!!$        rtask%bisVectorScalar     = .true.
!!$        rtask%ctask               =  PROBACTION_DUPLICATE
!!$        rtask%ssectionName        =  ssectionName
!!$        rtask%p_rvectorScalarSrc  => p_rproblemLevel%RvectorScalar(ivector)
!!$        rtask%p_rvectorScalarDest => rvector
!!$        rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
!!$        rtask%p_rproblemDest      => rproblemLevel%p_rproblem
!!$        rtask%p_rproblemLevelSrc  => p_rproblemLevel
!!$        rtask%p_rproblemLevelDest => rproblemLevel
!!$
!!$      else
!!$
!!$        !--- Block vector case -------------------------------------------------
!!$
!!$        ! Find problem name from which to duplicate entry of block vector (if any)
!!$        call parlst_getvalue_string(rparlist, ssectionName, 'vectorblock', sparameter, '')
!!$        call sys_toupper(sparameter)
!!$      
!!$        if (trim(sparameter) .ne. '') then
!!$          
!!$          ! Create new task for this scalar vector
!!$          nullify(rtask)
!!$          allocate(rtask)
!!$          nullify(rtask%p_rnextTask)
!!$          
!!$          ! Initialise optional parameters
!!$          sproblem = trim(rproblemLevel%p_rproblem%cproblem)
!!$          ilev     = rproblemLevel%ilev
!!$          ivector  = 1
!!$          iblock   = 1
!!$          
!!$          ! Get optional parameters if available
!!$          call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
!!$          
!!$          do itoken = 1,ntoken
!!$            call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
!!$            
!!$            iistart = 1
!!$            call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
!!$            call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
!!$            
!!$            if (trim(skeyword) .eq. 'PROBLEM') then
!!$              sproblem = trim(svalue)
!!$            elseif (trim(skeyword) .eq. 'ILEV') then
!!$              read(svalue,'(I10)') ilev
!!$            elseif (trim(skeyword) .eq. 'IVECTOR') then
!!$              read(svalue,'(I10)') ivector
!!$            elseif (trim(skeyword) .eq. 'IBLOCK') then
!!$              read(svalue,'(I10)') iblock
!!$            else
!!$              call output_line('Unsupported keyword: '//trim(skeyword),&
!!$                  OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$              call sys_halt()
!!$            end if
!!$          end do
!!$          
!!$          ! Find problem level in global problem structure
!!$          p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
!!$          if (.not.associated(p_rproblemLevel)) then
!!$            call output_line('Unable to find problem level in problem '//trim(sproblem),&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$            call sys_halt()
!!$          end if
!!$          
!!$          ! Get section name of block vector
!!$          call parlst_getvalue_string(rparlist, trim(sproblem),&
!!$              'vectorblock', ssubsectionName, isubstring=ivector)
!!$
!!$          ! Proceed with possibly prerequisite tasks
!!$          call problem_updateVecBl(p_rproblemLevel,&
!!$              p_rproblemLevel%RvectorBlock(ivector), rparlist,&
!!$              ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
!!$
!!$          ! Check if scalar subvector is available
!!$          p_rvectorBlock => p_rproblemLevel%RvectorBlock(ivector)
!!$          if ((size(p_rvectorBlock%RvectorBlock,1) .lt. iblock)) then
!!$            call output_line('Scalar subvector of block vector is not available',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$            call sys_halt()
!!$          end if
!!$
!!$          ! Specify new task for this cubature info structure
!!$          rtask%bisVectorScalar     = .true.
!!$          rtask%ctask               =  PROBACTION_DUPLICATE
!!$          rtask%ssectionName        =  ssectionName
!!$          rtask%p_rvectorScalarSrc  => p_rvectorBlock%RvectorBlock(iblock)
!!$          rtask%p_rvectorScalarDest => rvector
!!$          rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
!!$          rtask%p_rproblemDest      => rproblemLevel%p_rproblem
!!$          rtask%p_rproblemLevelSrc  => p_rproblemLevel
!!$          rtask%p_rproblemLevelDest => rproblemLevel
!!$
!!$        else
!!$          call output_line('Either vectorscalar of vectorblock must be present',&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateGroupFEMBlock')
!!$          call sys_halt()
!!$        end if
!!$      end if
!!$
!!$      ! Append task to task list (if not already present)
!!$      if (associated(p_rtasklist)) then
!!$        if (appendToTaskList(p_rtasklist, rtask)) then
!!$          deallocate(rtask)
!!$
!!$          ! That`s it we do not have to create this scalar vector
!!$          return
!!$        end if
!!$      else
!!$        p_rtasklist => rtask
!!$      end if
!!$
!!$      ! When should we perform this task?
!!$      call sys_countTokens(sperform, ntoken, ',', .false.)
!!$      istart = 1
!!$      do itoken=1,max(ntoken,1)
!!$        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
!!$        if (trim(adjustl(sperf)) .eq. 'INIT') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
!!$        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
!!$        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
!!$        else
!!$          call output_line('Unsupported perform type: '//trim(sperf),&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$      end do
!!$          
!!$      ! Get optional parameters
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'sdupstructure', sparameter, 'SHARE')
!!$      call sys_toupper(sparameter)
!!$      rtask%cdupStructure = get_dupType(sparameter)
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'sdupcontent', sparameter, 'SHARE')
!!$      call sys_toupper(sparameter)
!!$      rtask%cdupContent = get_dupType(sparameter)
!!$
!!$      ! Duplicate scalar vector
!!$      if (iand(rtask%iperform, iperform) .ne. 0) then
!!$        call lsyssc_duplicateVector(rtask%p_rvectorScalarSrc,&
!!$            rvector, rtask%cdupStructure, rtask%cdupContent)
!!$      end if
!!$      
    else
      call output_line('Unsupported action: '//trim(saction),&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateGroupFEMBlock')
      call sys_halt()
    end if
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Second part: generate content
!!$    !---------------------------------------------------------------------------
!!$
!!$    ! What creation method is adopted?
!!$    call parlst_getvalue_string(rparlist, ssectionName,&
!!$        'smethod', sparameter, '')
!!$    call sys_toupper(sparameter)
!!$    
!!$    if (trim(sparameter) .eq. 'EMPTY') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create empty vector
!!$      
!!$      call lsyssc_clearVector(rvector, 0.0_DP)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_EMPTY
!!$      
!!$    elseif (trim(sparameter) .eq. 'UNITY') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create unit vector
!!$      
!!$      call lsyssc_clearVector(rvector, 1.0_DP)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_UNITY
!!$      
!!$    elseif (trim(sparameter) .eq. 'LINF') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create vector by linearform evaluation         
!!$      
!!$      call assembleScalarVector(rparlist, rtask)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_LINF
!!$      
!!$    elseif (trim(sparameter) .eq. 'DOF') then
!!$
!!$      !-------------------------------------------------------------------------
!!$      ! Create coordinate vector for degrees of freedom
!!$      
!!$      call lin_calcDofCoords(rtask%p_rspatialDiscr, rvector)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_DOF
!!$
!!$    elseif (trim(sparameter) .eq. '') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create virtual vector, aka, do nothing
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_VIRTUAL
!!$      
!!$    else
!!$      call output_line('Unsupported method: '//trim(sparameter),&
!!$          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$      call sys_halt()
!!$    end if
!!$    
  contains

    ! Here, some auxiliary routines follow

    !***************************************************************************
    ! Search for given task in the list of tasks. If the task is not
    ! present in the list, then it is appended.
    function appendToTaskList(rtasklist, rtask) result(bexists)

      type(t_groupFEMBlockTask), intent(inout), target :: rtasklist
      type(t_groupFEMBlockTask), intent(in), target :: rtask
      
      ! local variable
      type(t_groupFEMBlockTask), pointer :: p_rtask
      logical :: bexists      

      p_rtask => rtasklist
      do while(associated(p_rtask))
        if (associated(p_rtask%p_rgroupFEMBlockDest,&
            rtask%p_rgroupFEMBlockDest)) then
          bexists = .true.
          return
        end if
        if (.not.associated(p_rtask%p_rnextTask)) then
          p_rtask%p_rnextTask => rtask
          bexists = .false.
          return
        else
          p_rtask => p_rtask%p_rnextTask
        end if
      end do
      
    end function appendToTaskList

    !**************************************************************
    ! Return the numeric value of the data type
    function get_dataType(sdataType) result(cdataType)

      character(len=*), intent(in) :: sdataType
      integer :: cdataType

      if (sdataType .eq. 'QUAD') then
        cdataType = ST_QUAD
      elseif (sdataType .eq. 'DOUBLE') then
        cdataType = ST_DOUBLE
      elseif (sdataType .eq. 'SINGLE') then
        cdataType = ST_SINGLE
      else
        cdataType = ST_NOHANDLE
      end if

    end function get_dataType
!!$
!!$    !**************************************************************
!!$    ! Return the numeric value of the duplication flad
!!$    function get_dupType(sdupType) result(cdupType)
!!$
!!$      character(len=*), intent(in) :: sdupType
!!$      integer :: cdupType
!!$
!!$      if (sdupType .eq. 'IGNORE') then
!!$        cdupType = LSYSSC_DUP_IGNORE
!!$      elseif (sdupType .eq. 'REMOVE') then
!!$        cdupType = LSYSSC_DUP_REMOVE
!!$      elseif (sdupType .eq. 'DISMISS') then
!!$        cdupType = LSYSSC_DUP_DISMISS
!!$      elseif (sdupType .eq. 'SHARE') then
!!$        cdupType = LSYSSC_DUP_SHARE
!!$      elseif (sdupType .eq. 'COPY') then
!!$        cdupType = LSYSSC_DUP_COPY
!!$      elseif (sdupType .eq. 'COPYOVERWRITE') then
!!$        cdupType = LSYSSC_DUP_COPYOVERWRITE
!!$      elseif (sdupType .eq. 'ASIS') then
!!$        cdupType = LSYSSC_DUP_ASIS
!!$      elseif (sdupType .eq. 'EMPTY') then
!!$        cdupType = LSYSSC_DUP_EMPTY
!!$      elseif (sdupType .eq. 'TEMPLATE') then
!!$        cdupType = LSYSSC_DUP_TEMPLATE
!!$      else
!!$        cdupType = 0
!!$      end if
!!$
!!$    end function get_dupType
!!$    
!!$    !**************************************************************
!!$    ! Assemble the scalar vector from bilinear form evaluation
!!$    subroutine assembleScalarVector(rparlist, rtask)
!!$
!!$      type(t_parlist), intent(in) :: rparlist
!!$      type(t_vectorTask), intent(inout) :: rtask
!!$
!!$      ! local variables
!!$      type(t_collection) :: rcollection
!!$      type(t_scalarCubatureInfo), pointer :: p_rcubInfo
!!$      character(len=PARLST_MLDATA) :: sparameter
!!$      integer :: i,j,itermCount
!!$
!!$      ! Get cubature info structure
!!$      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$          'cubatureinfo', sparameter)
!!$      call sys_toupper(sparameter)
!!$      p_rcubInfo => problem_getCubInfo(rtask%p_rproblemDest,&
!!$          rtask%p_rproblemLevelDest, sparameter)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%p_rcubatureInfo => p_rcubInfo
!!$      allocate(rtask%p_rform)
!!$      
!!$      ! Get number of terms in bilinear form
!!$      rtask%p_rform%itermCount = parlst_querysubstrings(&
!!$          rparlist, rtask%ssectionName, 'function')
!!$      
!!$      ! Get test and trial functions in bilinear form
!!$      do i = 1,rtask%p_rform%itermCount
!!$        call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$            'function', sparameter, isubString=i)
!!$        call sys_toupper(sparameter)
!!$        rtask%p_rform%Idescriptors(i) = der_igetID(sparameter)
!!$      end do
!!$      
!!$      ! Get constant(?) coefficients
!!$      itermCount = parlst_querysubstrings(rparlist,&
!!$          rtask%ssectionName, 'scoefficient')
!!$      
!!$      if (itermCount .gt. 0) then
!!$        itermCount = min(itermCount,rtask%p_rform%itermCount)
!!$        
!!$        ! Check if all coefficients are constant
!!$        do i=1,itermCount
!!$          call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$              'scoefficient', sparameter, isubString=i)
!!$          if (sys_isNumeric(sparameter)) then
!!$            call parlst_getvalue_double(rparlist, rtask%ssectionName,&
!!$                'scoefficient', rtask%p_rform%Dcoefficients(i), iarrayindex=i)
!!$          else
!!$            ! Not all coefficients are constant. Thus, we create a
!!$            ! function parser object and perform bilinear form
!!$            ! evaluation based on a generic callback function
!!$            allocate (rtask%p_rfparser)
!!$            call fparser_create(rtask%p_rfparser, itermCount)
!!$            
!!$            ! Fill function parser with data
!!$            do j=1,itermCount
!!$              call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$                  'scoefficient', sparameter, isubString=j)
!!$              call fparser_parseFunction(rtask%p_rfparser, j, trim(sparameter), (/'x','y','z'/))
!!$            end do
!!$            
!!$            ! That`s it
!!$            exit
!!$          end if
!!$        end do
!!$      else
!!$        rtask%p_rform%Dcoefficients(1:rtask%p_rform%itermCount) = 1.0_DP
!!$      end if
!!$
!!$      ! Assemble vector data from linear form
!!$      rcollection%p_rfparserQuickAccess1 => rtask%p_rfparser
!!$      if (associated(p_rcubInfo)) then
!!$        call linf_buildVectorScalar(rtask%p_rform, .true.,&
!!$            rtask%p_rvectorScalarDest, p_rcubInfo,&
!!$            problem_buildVectorSc_fparser_sim, rcollection)
!!$      elseif (associated(rtask%p_rspatialDiscr)) then
!!$        call linf_buildVectorScalar(rtask%p_rspatialDiscr,&
!!$            rtask%p_rform, .true., rtask%p_rvectorScalarDest,&
!!$            problem_buildVectorSc_fparser_sim, rcollection)
!!$      else
!!$        call output_line('Neither cubature info structure nor discretisation available',&
!!$            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$        call sys_halt()
!!$      end if
!!$      
!!$    end subroutine assembleScalarVector

  end subroutine problem_updateGroupFEMBlock1

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_updateGroupFEMBlock2(rtasklist, iperformSpec)

!<description>
    ! This subroutine updates the block of group finite element sets
    ! with the internal values assigned to the items of the task list.
!</description>

!<input>
    ! task list
    type(t_groupFEMBlockTask), intent(in), target :: rtasklist

    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>
!</subroutine>

!!$    ! local variables
!!$    type(t_collection) :: rcollection
!!$    type(t_vectorBlock) :: rvectorTemp
!!$    type(t_vectorTask), pointer :: p_rtask
!!$    integer :: iperform
!!$
!!$    iperform = PROBACTION_PERFORM_ALWAYS
!!$    if (present(iperformSpec)) iperform = iperformSpec
!!$    
!!$    ! Iterate over all tasks
!!$    p_rtask => rtasklist
!!$    taskloop: do while (associated(p_rtask))
!!$
!!$      ! Do we have to perform this task?
!!$      if (iand(p_rtask%iperform, iperform) .eq. 0) then
!!$        p_rtask => p_rtask%p_rnextTask
!!$        cycle taskloop
!!$      end if
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! First part: create or duplicate vector
!!$      !-------------------------------------------------------------------------
!!$      select case(p_rtask%ctask)
!!$      case(PROBACTION_CREATE)
!!$
!!$        ! Do we have to assemble a scalar vector?
!!$        if (p_rtask%bisVectorScalar) then
!!$          
!!$          !---------------------------------------------------------------------
!!$          ! Create scalar vector
!!$          
!!$          ! Do we have spatial discretisations?
!!$          if (.not.associated(p_rtask%p_rspatialDiscr)) then
!!$            call output_line('Unable to find spatial discretisation',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$            call sys_halt()
!!$          end if
!!$          
!!$          ! Clear preexisting vector
!!$          call lsyssc_releaseVector(p_rtask%p_rvectorScalarDest)
!!$          
!!$          ! Create vector structure and set data type
!!$          call lsyssc_createVector(p_rtask%p_rspatialDiscr,&
!!$              p_rtask%p_rvectorScalarDest, p_rtask%nvar, .false., p_rtask%cdataType)
!!$
!!$        else
!!$
!!$          !---------------------------------------------------------------------
!!$          ! Create block vector
!!$
!!$          ! Do we have block discretisations?
!!$          if (.not.associated(p_rtask%p_rblockDiscr)) then
!!$            call output_line('Unable to find block discretisation',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$            call sys_halt()
!!$          end if
!!$
!!$          ! Clear preexisting block vector
!!$          call lsysbl_releaseVector(p_rtask%p_rvectorBlockDest)
!!$
!!$          ! Create block vector from discretisation
!!$          call lsysbl_createVecBlockByDiscr(p_rtask%p_rblockDiscr,&
!!$              p_rtask%p_rvectorBlockDest)
!!$          
!!$          ! Proceed with next task
!!$          p_rtask => p_rtask%p_rnextTask
!!$          
!!$          ! That`s it for block vectors
!!$          cycle taskloop
!!$        end if
!!$
!!$      case(PROBACTION_DUPLICATE)
!!$
!!$        !-----------------------------------------------------------------------
!!$        ! Duplicate scalar/block vector
!!$
!!$        if (p_rtask%bisVectorScalar) then
!!$          ! Duplicate scalar vector
!!$          call lsyssc_duplicateVector(p_rtask%p_rvectorScalarSrc,&
!!$              p_rtask%p_rvectorScalarDest, p_rtask%cdupStructure, p_rtask%cdupContent)
!!$        else
!!$          if (associated(p_rtask%p_rvectorScalarSrc)) then
!!$            ! Duplicate scalar vector as 1-block vector
!!$            call lsysbl_createVecFromScalar(p_rtask%p_rvectorScalarSrc, rvectorTemp)
!!$
!!$            ! Copy content of temporal vector by hand
!!$            p_rtask%p_rvectorBlockDest%NEQ         = rvectorTemp%NEQ
!!$            p_rtask%p_rvectorBlockDest%nblocks     = 1
!!$
!!$            allocate(p_rtask%p_rvectorBlockDest%RvectorBlock(1))
!!$            p_rtask%p_rvectorBlockDest%RvectorBlock(1) = rvectorTemp%RvectorBlock(1)
!!$
!!$            ! Release temporal vector
!!$            call lsysbl_releaseVector(rvectorTemp)
!!$
!!$          elseif(associated(p_rtask%p_rvectorBlockSrc)) then
!!$            ! Duplicate block vector
!!$            call lsysbl_duplicateVector(p_rtask%p_rvectorBlockSrc,&
!!$                p_rtask%p_rvectorBlockDest, p_rtask%cdupStructure, p_rtask%cdupContent)
!!$          else
!!$            call output_line('Neither scalar nor block vector can be duplicated',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$            call sys_halt()
!!$          end if
!!$        end if
!!$
!!$      case default
!!$        call output_line('Unsupported action: '//sys_siL(p_rtask%ctask,3),&
!!$            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMamtrix')
!!$        call sys_halt()
!!$      end select
!!$
!!$      !-------------------------------------------------------------------------
!!$      ! Second part: generate content (only for scalar vectors)
!!$      !-------------------------------------------------------------------------
!!$
!!$      ! What creation method is adopted?
!!$      select case (p_rtask%cmethod)
!!$      case(PROBMETH_VECTOR_EMPTY)
!!$
!!$        !---------------------------------------------------------------------
!!$        ! Create empty vector
!!$
!!$        call lsyssc_clearVector(p_rtask%p_rvectorScalarDest, 0.0_DP)
!!$
!!$      case(PROBMETH_VECTOR_UNITY)
!!$
!!$        !---------------------------------------------------------------------
!!$        ! Create unit vector
!!$
!!$        call lsyssc_clearVector(p_rtask%p_rvectorScalarDest, 1.0_DP)
!!$                
!!$      case(PROBMETH_VECTOR_LINF)
!!$
!!$        !---------------------------------------------------------------------
!!$        ! Create vector by bilinearform evaluation         
!!$
!!$        rcollection%p_rfparserQuickAccess1 => p_rtask%p_rfparser
!!$        if (associated(p_rtask%p_rcubatureInfo)) then
!!$          call linf_buildVectorScalar(p_rtask%p_rform, .true.,&
!!$              p_rtask%p_rvectorScalarDest, p_rtask%p_rcubatureInfo,&
!!$              problem_buildVectorSc_fparser_sim, rcollection)
!!$        elseif (associated(p_rtask%p_rspatialDiscr)) then
!!$          call linf_buildVectorScalar(p_rtask%p_rspatialDiscr,&
!!$              p_rtask%p_rform, .true., p_rtask%p_rvectorScalarDest,&
!!$              problem_buildVectorSc_fparser_sim, rcollection)
!!$        else
!!$          call output_line('Neither cubature info structure nor discretisation available',&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$          call sys_halt()
!!$        end if
!!$        
!!$      case(PROBMETH_VECTOR_VIRTUAL)
!!$        
!!$        !---------------------------------------------------------------------
!!$        ! Create virtual vector, aka, do nothing
!!$        
!!$      case default
!!$        call output_line('Unsupported method: '//sys_siL(p_rtask%cmethod,3),&
!!$            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$        call sys_halt()
!!$      end select
!!$      
!!$      ! Proceed with next task
!!$      p_rtask => p_rtask%p_rnextTask
!!$    end do taskloop
    
  end subroutine problem_updateGroupFEMBlock2

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_clearGroupFEMBlockTL(p_rtasklist)

!<description>
    ! This subroutine clears the block of group finite element sets
    ! task list.
!</description>

!<input>
    ! task list
    type(t_groupFEMBlockTask), pointer :: p_rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_groupFEMBlockTask), pointer :: p_rtask,p_rtaskSave

    p_rtask => p_rtasklist
    do while(associated(p_rtask))
      p_rtaskSave => p_rtask
      p_rtask     => p_rtask%p_rnextTask
      deallocate(p_rtaskSave)
    end do
    
  end subroutine problem_clearGroupFEMBlockTL

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_infoGroupFEMBlockTL(rtasklist)

!<description>
    ! This subroutine outputs information about the task list of
    ! blocks of group finite element sets.
!</description>

!<input>
    ! task list
    type(t_groupFEMBlockTask), intent(in), target :: rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_groupFEMBlockTask), pointer :: p_rtask

    p_rtask => rtasklist
    do while(associated(p_rtask))

      select case(p_rtask%ctask)
      case (PROBACTION_NONE)
        call output_lbrk()
        call output_line ('GroupFEMBlockTask: DO NOTHING')

      case (PROBACTION_CREATE)
        call output_lbrk()
        call output_line ('GroupFEMBlockTask: CREATE')

      case (PROBACTION_DUPLICATE)
        call output_lbrk()
        call output_line ('GroupFEMBlockTask: DUPLICATE')

      case default
        call output_line ('GroupFEMBockTask: UNSUPPORTED')
      end select

!!$      call output_line ('-----------------------')
!!$      call output_line ('vector type               : '//merge('SCALAR', 'BLOCK ', p_rtask%bisVectorScalar))
!!$      call output_line ('section name              : '//trim(p_rtask%ssectionName))
!!$      call output_line ('destination problem       : '//trim(p_rtask%p_rproblemDest%cproblem))
!!$      call output_line ('destination problem level : '//trim(sys_siL(p_rtask%p_rproblemLevelDest%ilev,15)))
!!$      call output_line ('destination vector scalar : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorScalarDest)))
!!$      call output_line ('destination vector block  : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorBlockDest)))
!!$      if (associated(p_rtask%p_rproblemSrc))&
!!$      call output_line ('source problem            : '//trim(p_rtask%p_rproblemSrc%cproblem))
!!$      if (associated(p_rtask%p_rproblemLevelSrc))&
!!$      call output_line ('source problem level      : '//trim(sys_siL(p_rtask%p_rproblemLevelSrc%ilev,15)))
!!$      call output_line ('source vector scalar      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorScalarSrc)))
!!$      call output_line ('source vector block       : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorBlockSrc)))
!!$      call output_line ('iperform                  : '//trim(sys_siL(int(p_rtask%iperform),15)))
!!$      call output_line ('method                    : '//trim(sys_siL(p_rtask%cmethod,15)))
!!$      call output_line ('data type                 : '//trim(sys_siL(p_rtask%cdatatype,15)))
!!$      call output_line ('nvar                      : '//trim(sys_siL(p_rtask%nvar,15)))
!!$      call output_line ('spatial discretisation:     '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rspatialDiscr)))
!!$      call output_line ('block discretisation      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rblockDiscr)))
!!$      call output_line ('cubature info structure   : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rcubatureInfo)))
!!$      call output_line ('function parser           : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rfparser)))
!!$      call output_line ('linear form evaluator     : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rform)))

      p_rtask => p_rtask%p_rnextTask
    end do
    
  end subroutine problem_infoGroupFEMBlockTL

  !*****************************************************************************

!<subroutine>

  subroutine problem_initAFCstabAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises all stabilisation structures of
    ! AFC-type of the given problem level with the values provided by
    ! the parameter list.  If the optional parameter p_rtasklist is
    ! given, then the task list which is created during the
    ! initialisation is returned.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! OPTIONAL: task list
    type(t_afcstabTask), intent(inout), pointer, optional :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    type(t_afcstabTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update all stabilisation structures of AFC-type
    call problem_updateAFCstabAll(rproblemLevel,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearAFCstabTL(p_rtasklistLocal)
    end if

  end subroutine problem_initAFCstabAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_updateAFCstabAll(rproblemLevel,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine the task list for all stabilisation structures
    ! of AFC-type of the given problem level with the values provided
    ! by the parameter list.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! pointer to the task list
    type(t_afcstabTask), pointer :: p_rtasklist
!</intputoutput>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: ssubsectionName
    integer :: i,nafcstab

    ! Consistency check
    nafcstab = parlst_querysubstrings(rparlist, ssectionName,&
                                      'afcstab')
    if ((nafcstab .gt. 0) .and.&
        (nafcstab .ne. size(rproblemLevel%Rafcstab))) then
      call output_line('Invalid number of stabilisation structures of AFC-type',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateAFCstabAll')
      call sys_halt()
    end if
   
    ! Update all stabilisation structures of AFC-type
    do i=1,nafcstab
      call parlst_getvalue_string(rparlist, ssectionName,&
          'afcstab', ssubsectionName, isubstring=i)
      call problem_updateAFCstab1(rproblemLevel,&
          rproblemLevel%Rafcstab(i), rparlist, ssubsectionName,&
          p_rtasklist, rproblemTopLevel, iperformSpec)
    end do
    
  end subroutine problem_updateAFCstabAll

  !*****************************************************************************

!<subroutine>

  subroutine problem_initAFCstab(rproblemLevel, rafcstab,&
      rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine initialises the given stabilisation structure of
    ! AFC-type on the given problem level with the values provided by
    ! the parameter list. If the optional parameter p_rtasklist is
    ! given, then the task list which is created during the
    ! initialisation is returned.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which stabilisation
    ! structures of AFC-type need to be created and copied
    type(t_afcstabTask), pointer, optional :: p_rtasklist
!</inputoutput>

!<output>
    ! stabilisation structure of AFC-type
    type(t_afcstab), intent(out), target :: rafcstab
!</output>
!</subroutine>

    ! local variables
    type(t_afcstabTask), pointer :: p_rtasklistLocal => null()

    if (present(p_rtasklist)) p_rtasklistLocal => p_rtasklist
    
    ! Update given stabilisation structure of AFC-type
    call problem_updateAFCstab1(rproblemLevel, rafcstab,&
        rparlist, ssectionName, p_rtasklistLocal, rproblemTopLevel, iperformSpec)

    ! Release temporal task list if required
    if (.not.present(p_rtasklist)) then
      call problem_clearAFCstabTL(p_rtasklistLocal)
    end if

  end subroutine problem_initAFCstab

  !*****************************************************************************

!<subroutine>

  recursive subroutine problem_updateAFCstab1(rproblemLevel,&
      rafcstab, rparlist, ssectionName, p_rtasklist, rproblemTopLevel, iperformSpec)

!<description>
    ! This subroutine generates the task list to initialise/update the
    ! given stabilisation structure of AFC-type on the given problem
    ! level with the values provided by the parameter list.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: top-level problem structure
    ! If not present, then the problem associated with the problem
    ! level structure serves as top-level problem
    type(t_problem), intent(in), target, optional :: rproblemTopLevel

    ! OPTIONAL: specifier for the tasks to be performed
    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>

!<inputoutput>
    ! task list which represents the order in which stabilisation
    ! structures of AFC-type need to be created and copied
    type(t_afcstabTask), pointer :: p_rtasklist

    ! stabilisation structure of AFC-type
    type(t_afcstab), intent(inout), target :: rafcstab
!</inputoutput>
!</subroutine>

!!$    ! local variables
!!$    type(t_vectorBlock), pointer :: p_rvectorBlock
!!$    type(t_vectorTask), pointer :: rtask
!!$    type(t_problem), pointer :: p_rproblemTopLevel
!!$    type(t_problemLevel), pointer :: p_rproblemLevel
!!$    character(len=PARLST_MLDATA) :: skeyword,sparameter,sproblem,stoken,svalue
!!$    character(len=SYS_STRLEN) :: ssubsectionName
!!$    character(len=SYS_STRLEN) :: saction,sperform,sperf
!!$    integer :: iblock,iistart,ilev,ivector,istart,itoken,ntoken
!!$    integer(i32) :: iperform
!!$
!!$    ! Set pointer to problem structure
!!$    p_rproblemTopLevel => rproblemLevel%p_rproblem
!!$    if (present(rproblemTopLevel)) p_rproblemTopLevel => rproblemTopLevel
!!$
!!$    iperform = PROBACTION_PERFORM_ALWAYS
!!$    if (present(iperformSpec)) iperform = iperformSpec
!!$
!!$    ! Get type of action
!!$    call parlst_getvalue_string(rparlist, ssectionName, 'saction', saction)
!!$    call parlst_getvalue_string(rparlist, ssectionName, 'sperform', sperform)
!!$    call sys_toupper(saction)
!!$    call sys_toupper(sperform)
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! First part: generate structure
!!$    !---------------------------------------------------------------------------
!!$
!!$    ! What action should be performed?
!!$    if (trim(saction) .eq. 'NONE') then
!!$      !-------------------------------------------------------------------------
!!$      ! Do nothing
!!$
!!$    elseif (trim(saction) .eq. 'CREATE') then
!!$
!!$      !-------------------------------------------------------------------------
!!$      ! Create scalar vector
!!$      !
!!$      ! SYNTAX: vectorscalar(n) =
!!$      !           VectorScalar1
!!$      !                ...
!!$      !           VectorScalar2
!!$
!!$      ! Create new task
!!$      nullify(rtask)
!!$      allocate(rtask)
!!$      nullify(rtask%p_rnextTask)
!!$      rtask%bisVectorScalar     = .true.
!!$      rtask%ctask               =  PROBACTION_CREATE
!!$      rtask%ssectionName        =  ssectionName
!!$      rtask%p_rvectorScalarDest => rvector
!!$      rtask%p_rproblemDest      => rproblemLevel%p_rproblem
!!$      rtask%p_rproblemLevelDest => rproblemLevel
!!$
!!$      ! Append task to task list (if not already present)
!!$      if (associated(p_rtasklist)) then
!!$        if (appendToTaskList(p_rtasklist, rtask)) then
!!$          deallocate(rtask)
!!$
!!$          ! That`s it we do not have to create this scalar vector
!!$          return
!!$        end if
!!$      else
!!$        p_rtasklist => rtask
!!$      end if
!!$
!!$      ! When should we perform this task?
!!$      call sys_countTokens(sperform, ntoken, ',', .false.)
!!$      istart = 1
!!$      do itoken=1,max(ntoken,1)
!!$        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
!!$        if (trim(adjustl(sperf)) .eq. 'INIT') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
!!$        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
!!$        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
!!$        else
!!$          call output_line('Unsupported perform type: '//trim(sperf),&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$      end do
!!$      
!!$      ! Get optional parameters
!!$    call parlst_getvalue_int(rparlist, ssectionName,&
!!$        'cdatatype', rtask%cdataType, ST_NOHANDLE)
!!$    if (rtask%cdataType .eq. ST_NOHANDLE) then
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'sdatatype', sparameter, 'DOUBLE')
!!$      call sys_toupper(sparameter)
!!$      rtask%cdataType = get_dataType(trim(sparameter))
!!$    end if
!!$      call parlst_getvalue_int(rparlist, ssectionName,&
!!$          'nvar', rtask%nvar, 1)
!!$      
!!$      ! Get spatial discretisation for test and trial space
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'discretisation', sparameter)
!!$      call sys_toupper(sparameter)
!!$      rtask%p_rspatialDiscr => problem_getSpatialDiscr(&
!!$          rproblemLevel%p_rproblem, rproblemLevel, sparameter)
!!$      
!!$      ! Should we perform this task now?
!!$      if (iand(rtask%iperform, iperform) .ne. 0) then
!!$        ! Do we have spatial discretisations?
!!$        if (.not.associated(rtask%p_rspatialDiscr)) then
!!$          call output_line('Unable to find spatial discretisation',&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$        
!!$        ! Clear preexisting scalar vector
!!$        call lsyssc_releaseVector(rvector)
!!$        
!!$        ! Create vector structure and set data type
!!$        call lsyssc_createVector(rtask%p_rspatialDiscr, rvector,&
!!$            rtask%nvar, .false., rtask%cdataType)
!!$      end if
!!$
!!$    elseif (trim(saction) .eq. 'DUPLICATE') then
!!$
!!$      !-------------------------------------------------------------------------
!!$      ! Duplicate scalar vector or subvector of a block vector
!!$      !
!!$      ! SYNTAX: vectorscalar = name,ivector:#,ilev:#,...
!!$      !     OR  vectorblock  = name,ivector:#,ilev:#,iblock:#
!!$      !
!!$      ! If ilev is not given then the level of the current problem
!!$      ! level structure is adopted. If iblock are not given then
!!$      ! standard value 1 is used in both cases
!!$
!!$      ! Find problem name from which to duplicate scalar vector (if any)
!!$      call parlst_getvalue_string(rparlist, ssectionName, 'vectorscalar', sparameter, '')
!!$      call sys_toupper(sparameter)
!!$
!!$      if (trim(sparameter) .ne. '') then
!!$        
!!$        !--- Scalar vector case ------------------------------------------------
!!$
!!$        ! Create new task for this scalar vector
!!$        nullify(rtask)
!!$        allocate(rtask)
!!$        nullify(rtask%p_rnextTask)
!!$        
!!$        ! Initialise optional parameters
!!$        sproblem = trim(rproblemLevel%p_rproblem%cproblem)
!!$        ilev     = rproblemLevel%ilev
!!$        ivector  = 1
!!$        
!!$        ! Get optional parameters if available
!!$        call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
!!$        
!!$        do itoken = 1,ntoken
!!$          call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
!!$          
!!$          iistart = 1
!!$          call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
!!$          call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
!!$          
!!$          if (trim(skeyword) .eq. 'PROBLEM') then
!!$            sproblem = trim(svalue)
!!$          elseif (trim(skeyword) .eq. 'ILEV') then
!!$            read(svalue,'(I10)') ilev
!!$          elseif (trim(skeyword) .eq. 'IVECTOR') then
!!$            read(svalue,'(I10)') ivector
!!$          else
!!$            call output_line('Unsupported keyword: '//trim(skeyword),&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$            call sys_halt()
!!$          end if
!!$        end do
!!$        
!!$        ! Find problem level in global problem structure
!!$        p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
!!$        if (.not.associated(p_rproblemLevel)) then
!!$          call output_line('Unable to find problem level in problem '//trim(sproblem),&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$        
!!$        ! Get section name of scalar vector
!!$        call parlst_getvalue_string(rparlist, trim(sproblem),&
!!$            'vectorscalar', ssubsectionName, isubstring=ivector)
!!$
!!$        ! Proceed with possibly prerequisite tasks
!!$        call problem_updateVectorScalar(p_rproblemLevel,&
!!$            p_rproblemLevel%RvectorScalar(ivector), rparlist,&
!!$            ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
!!$        
!!$        ! Specify new task for this cubature info structure
!!$        rtask%bisVectorScalar     = .true.
!!$        rtask%ctask               =  PROBACTION_DUPLICATE
!!$        rtask%ssectionName        =  ssectionName
!!$        rtask%p_rvectorScalarSrc  => p_rproblemLevel%RvectorScalar(ivector)
!!$        rtask%p_rvectorScalarDest => rvector
!!$        rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
!!$        rtask%p_rproblemDest      => rproblemLevel%p_rproblem
!!$        rtask%p_rproblemLevelSrc  => p_rproblemLevel
!!$        rtask%p_rproblemLevelDest => rproblemLevel
!!$
!!$      else
!!$
!!$        !--- Block vector case -------------------------------------------------
!!$
!!$        ! Find problem name from which to duplicate entry of block vector (if any)
!!$        call parlst_getvalue_string(rparlist, ssectionName, 'vectorblock', sparameter, '')
!!$        call sys_toupper(sparameter)
!!$      
!!$        if (trim(sparameter) .ne. '') then
!!$          
!!$          ! Create new task for this scalar vector
!!$          nullify(rtask)
!!$          allocate(rtask)
!!$          nullify(rtask%p_rnextTask)
!!$          
!!$          ! Initialise optional parameters
!!$          sproblem = trim(rproblemLevel%p_rproblem%cproblem)
!!$          ilev     = rproblemLevel%ilev
!!$          ivector  = 1
!!$          iblock   = 1
!!$          
!!$          ! Get optional parameters if available
!!$          call sys_countTokens(sparameter, ntoken, ',', .false.); istart=1
!!$          
!!$          do itoken = 1,ntoken
!!$            call sys_getNextToken (sparameter, stoken, istart, ',', .false.)
!!$            
!!$            iistart = 1
!!$            call sys_getNextToken (stoken, skeyword, iistart, ":", .false.)
!!$            call sys_getNextToken (stoken, svalue, iistart, ":", .false.)
!!$            
!!$            if (trim(skeyword) .eq. 'PROBLEM') then
!!$              sproblem = trim(svalue)
!!$            elseif (trim(skeyword) .eq. 'ILEV') then
!!$              read(svalue,'(I10)') ilev
!!$            elseif (trim(skeyword) .eq. 'IVECTOR') then
!!$              read(svalue,'(I10)') ivector
!!$            elseif (trim(skeyword) .eq. 'IBLOCK') then
!!$              read(svalue,'(I10)') iblock
!!$            else
!!$              call output_line('Unsupported keyword: '//trim(skeyword),&
!!$                  OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$              call sys_halt()
!!$            end if
!!$          end do
!!$          
!!$          ! Find problem level in global problem structure
!!$          p_rproblemLevel => problem_getLevel(p_rproblemTopLevel, trim(sproblem), ilev)
!!$          if (.not.associated(p_rproblemLevel)) then
!!$            call output_line('Unable to find problem level in problem '//trim(sproblem),&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$            call sys_halt()
!!$          end if
!!$          
!!$          ! Get section name of block vector
!!$          call parlst_getvalue_string(rparlist, trim(sproblem),&
!!$              'vectorblock', ssubsectionName, isubstring=ivector)
!!$
!!$          ! Proceed with possibly prerequisite tasks
!!$          call problem_updateVecBl(p_rproblemLevel,&
!!$              p_rproblemLevel%RvectorBlock(ivector), rparlist,&
!!$              ssubsectionName, p_rtasklist, p_rproblemTopLevel, iperformSpec)
!!$
!!$          ! Check if scalar subvector is available
!!$          p_rvectorBlock => p_rproblemLevel%RvectorBlock(ivector)
!!$          if ((size(p_rvectorBlock%RvectorBlock,1) .lt. iblock)) then
!!$            call output_line('Scalar subvector of block vector is not available',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$            call sys_halt()
!!$          end if
!!$
!!$          ! Specify new task for this cubature info structure
!!$          rtask%bisVectorScalar     = .true.
!!$          rtask%ctask               =  PROBACTION_DUPLICATE
!!$          rtask%ssectionName        =  ssectionName
!!$          rtask%p_rvectorScalarSrc  => p_rvectorBlock%RvectorBlock(iblock)
!!$          rtask%p_rvectorScalarDest => rvector
!!$          rtask%p_rproblemSrc       => p_rproblemLevel%p_rproblem
!!$          rtask%p_rproblemDest      => rproblemLevel%p_rproblem
!!$          rtask%p_rproblemLevelSrc  => p_rproblemLevel
!!$          rtask%p_rproblemLevelDest => rproblemLevel
!!$
!!$        else
!!$          call output_line('Either vectorscalar of vectorblock must be present',&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$      end if
!!$
!!$      ! Append task to task list (if not already present)
!!$      if (associated(p_rtasklist)) then
!!$        if (appendToTaskList(p_rtasklist, rtask)) then
!!$          deallocate(rtask)
!!$
!!$          ! That`s it we do not have to create this scalar vector
!!$          return
!!$        end if
!!$      else
!!$        p_rtasklist => rtask
!!$      end if
!!$
!!$      ! When should we perform this task?
!!$      call sys_countTokens(sperform, ntoken, ',', .false.)
!!$      istart = 1
!!$      do itoken=1,max(ntoken,1)
!!$        call sys_getNextToken (sperform, sperf, istart, ',', .false.)
!!$        if (trim(adjustl(sperf)) .eq. 'INIT') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_INIT)
!!$        elseif (trim(adjustl(sperf)) .eq. 'UPDATE') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_UPDATE)
!!$        elseif (trim(adjustl(sperf)) .eq. 'ALWAYS') then
!!$          rtask%iperform = ior(rtask%iperform, PROBACTION_PERFORM_ALWAYS)
!!$        else
!!$          call output_line('Unsupported perform type: '//trim(sperf),&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$          call sys_halt()
!!$        end if
!!$      end do
!!$          
!!$      ! Get optional parameters
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'sdupstructure', sparameter, 'SHARE')
!!$      call sys_toupper(sparameter)
!!$      rtask%cdupStructure = get_dupType(sparameter)
!!$      call parlst_getvalue_string(rparlist, ssectionName,&
!!$          'sdupcontent', sparameter, 'SHARE')
!!$      call sys_toupper(sparameter)
!!$      rtask%cdupContent = get_dupType(sparameter)
!!$
!!$      ! Duplicate scalar vector
!!$      if (iand(rtask%iperform, iperform) .ne. 0) then
!!$        call lsyssc_duplicateVector(rtask%p_rvectorScalarSrc,&
!!$            rvector, rtask%cdupStructure, rtask%cdupContent)
!!$      end if
!!$      
!!$    else
!!$      call output_line('Unsupported action: '//trim(saction),&
!!$          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$      call sys_halt()
!!$    end if
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Second part: generate content
!!$    !---------------------------------------------------------------------------
!!$
!!$    ! What creation method is adopted?
!!$    call parlst_getvalue_string(rparlist, ssectionName,&
!!$        'smethod', sparameter, '')
!!$    call sys_toupper(sparameter)
!!$    
!!$    if (trim(sparameter) .eq. 'EMPTY') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create empty vector
!!$      
!!$      call lsyssc_clearVector(rvector, 0.0_DP)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_EMPTY
!!$      
!!$    elseif (trim(sparameter) .eq. 'UNITY') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create unit vector
!!$      
!!$      call lsyssc_clearVector(rvector, 1.0_DP)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_UNITY
!!$      
!!$    elseif (trim(sparameter) .eq. 'LINF') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create vector by linearform evaluation         
!!$      
!!$      call assembleScalarVector(rparlist, rtask)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_LINF
!!$      
!!$    elseif (trim(sparameter) .eq. 'DOF') then
!!$
!!$      !-------------------------------------------------------------------------
!!$      ! Create coordinate vector for degrees of freedom
!!$      
!!$      call lin_calcDofCoords(rtask%p_rspatialDiscr, rvector)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_DOF
!!$
!!$    elseif (trim(sparameter) .eq. '') then
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! Create virtual vector, aka, do nothing
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%cmethod = PROBMETH_VECTOR_VIRTUAL
!!$      
!!$    else
!!$      call output_line('Unsupported method: '//trim(sparameter),&
!!$          OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$      call sys_halt()
!!$    end if
!!$    
!!$  contains
!!$
!!$    ! Here, some auxiliary routines follow
!!$
!!$    !***************************************************************************
!!$    ! Search for given task in the list of tasks. If the task is not
!!$    ! present in the list, then it is appended.
!!$    function appendToTaskList(rtasklist, rtask) result(bexists)
!!$
!!$      type(t_vectorTask), intent(inout), target :: rtasklist
!!$      type(t_vectorTask), intent(in), target :: rtask
!!$      
!!$      ! local variable
!!$      type(t_vectorTask), pointer :: p_rtask
!!$      logical :: bexists      
!!$
!!$      p_rtask => rtasklist
!!$      do while(associated(p_rtask))
!!$        if (associated(p_rtask%p_rvectorScalarDest,&
!!$            rtask%p_rvectorScalarDest)) then
!!$          bexists = .true.
!!$          return
!!$        end if
!!$        if (.not.associated(p_rtask%p_rnextTask)) then
!!$          p_rtask%p_rnextTask => rtask
!!$          bexists = .false.
!!$          return
!!$        else
!!$          p_rtask => p_rtask%p_rnextTask
!!$        end if
!!$      end do
!!$      
!!$    end function appendToTaskList
!!$
!!$    !**************************************************************
!!$    ! Return the numeric value of the data type
!!$    function get_dataType(sdataType) result(cdataType)
!!$
!!$      character(len=*), intent(in) :: sdataType
!!$      integer :: cdataType
!!$
!!$      if (sdataType .eq. 'QUAD') then
!!$        cdataType = ST_QUAD
!!$      elseif (sdataType .eq. 'DOUBLE') then
!!$        cdataType = ST_DOUBLE
!!$      elseif (sdataType .eq. 'SINGLE') then
!!$        cdataType = ST_SINGLE
!!$      else
!!$        cdataType = ST_NOHANDLE
!!$      end if
!!$
!!$    end function get_dataType
!!$
!!$    !**************************************************************
!!$    ! Return the numeric value of the duplication flad
!!$    function get_dupType(sdupType) result(cdupType)
!!$
!!$      character(len=*), intent(in) :: sdupType
!!$      integer :: cdupType
!!$
!!$      if (sdupType .eq. 'IGNORE') then
!!$        cdupType = LSYSSC_DUP_IGNORE
!!$      elseif (sdupType .eq. 'REMOVE') then
!!$        cdupType = LSYSSC_DUP_REMOVE
!!$      elseif (sdupType .eq. 'DISMISS') then
!!$        cdupType = LSYSSC_DUP_DISMISS
!!$      elseif (sdupType .eq. 'SHARE') then
!!$        cdupType = LSYSSC_DUP_SHARE
!!$      elseif (sdupType .eq. 'COPY') then
!!$        cdupType = LSYSSC_DUP_COPY
!!$      elseif (sdupType .eq. 'COPYOVERWRITE') then
!!$        cdupType = LSYSSC_DUP_COPYOVERWRITE
!!$      elseif (sdupType .eq. 'ASIS') then
!!$        cdupType = LSYSSC_DUP_ASIS
!!$      elseif (sdupType .eq. 'EMPTY') then
!!$        cdupType = LSYSSC_DUP_EMPTY
!!$      elseif (sdupType .eq. 'TEMPLATE') then
!!$        cdupType = LSYSSC_DUP_TEMPLATE
!!$      else
!!$        cdupType = 0
!!$      end if
!!$
!!$    end function get_dupType
!!$    
!!$    !**************************************************************
!!$    ! Assemble the scalar vector from bilinear form evaluation
!!$    subroutine assembleScalarVector(rparlist, rtask)
!!$
!!$      type(t_parlist), intent(in) :: rparlist
!!$      type(t_vectorTask), intent(inout) :: rtask
!!$
!!$      ! local variables
!!$      type(t_collection) :: rcollection
!!$      type(t_scalarCubatureInfo), pointer :: p_rcubInfo
!!$      character(len=PARLST_MLDATA) :: sparameter
!!$      integer :: i,j,itermCount
!!$
!!$      ! Get cubature info structure
!!$      call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$          'cubatureinfo', sparameter)
!!$      call sys_toupper(sparameter)
!!$      p_rcubInfo => problem_getCubInfo(rtask%p_rproblemDest,&
!!$          rtask%p_rproblemLevelDest, sparameter)
!!$      
!!$      ! Update internal data of the task item
!!$      rtask%p_rcubatureInfo => p_rcubInfo
!!$      allocate(rtask%p_rform)
!!$      
!!$      ! Get number of terms in bilinear form
!!$      rtask%p_rform%itermCount = parlst_querysubstrings(&
!!$          rparlist, rtask%ssectionName, 'function')
!!$      
!!$      ! Get test and trial functions in bilinear form
!!$      do i = 1,rtask%p_rform%itermCount
!!$        call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$            'function', sparameter, isubString=i)
!!$        call sys_toupper(sparameter)
!!$        rtask%p_rform%Idescriptors(i) = der_igetID(sparameter)
!!$      end do
!!$      
!!$      ! Get constant(?) coefficients
!!$      itermCount = parlst_querysubstrings(rparlist,&
!!$          rtask%ssectionName, 'scoefficient')
!!$      
!!$      if (itermCount .gt. 0) then
!!$        itermCount = min(itermCount,rtask%p_rform%itermCount)
!!$        
!!$        ! Check if all coefficients are constant
!!$        do i=1,itermCount
!!$          call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$              'scoefficient', sparameter, isubString=i)
!!$          if (sys_isNumeric(sparameter)) then
!!$            call parlst_getvalue_double(rparlist, rtask%ssectionName,&
!!$                'scoefficient', rtask%p_rform%Dcoefficients(i), iarrayindex=i)
!!$          else
!!$            ! Not all coefficients are constant. Thus, we create a
!!$            ! function parser object and perform bilinear form
!!$            ! evaluation based on a generic callback function
!!$            allocate (rtask%p_rfparser)
!!$            call fparser_create(rtask%p_rfparser, itermCount)
!!$            
!!$            ! Fill function parser with data
!!$            do j=1,itermCount
!!$              call parlst_getvalue_string(rparlist, rtask%ssectionName,&
!!$                  'scoefficient', sparameter, isubString=j)
!!$              call fparser_parseFunction(rtask%p_rfparser, j, trim(sparameter), (/'x','y','z'/))
!!$            end do
!!$            
!!$            ! That`s it
!!$            exit
!!$          end if
!!$        end do
!!$      else
!!$        rtask%p_rform%Dcoefficients(1:rtask%p_rform%itermCount) = 1.0_DP
!!$      end if
!!$
!!$      ! Assemble vector data from linear form
!!$      rcollection%p_rfparserQuickAccess1 => rtask%p_rfparser
!!$      if (associated(p_rcubInfo)) then
!!$        call linf_buildVectorScalar(rtask%p_rform, .true.,&
!!$            rtask%p_rvectorScalarDest, p_rcubInfo,&
!!$            problem_buildVectorSc_fparser_sim, rcollection)
!!$      elseif (associated(rtask%p_rspatialDiscr)) then
!!$        call linf_buildVectorScalar(rtask%p_rspatialDiscr,&
!!$            rtask%p_rform, .true., rtask%p_rvectorScalarDest,&
!!$            problem_buildVectorSc_fparser_sim, rcollection)
!!$      else
!!$        call output_line('Neither cubature info structure nor discretisation available',&
!!$            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVectorScalar')
!!$        call sys_halt()
!!$      end if
!!$      
!!$    end subroutine assembleScalarVector

  end subroutine problem_updateAFCstab1

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_updateAFCstab2(rtasklist, iperformSpec)

!<description>
    ! This subroutine updates the stabilisation structure of AFC-type
    ! with the internal values assigned to the items of the task list.
!</description>

!<input>
    ! task list
    type(t_afcstabTask), intent(in), target :: rtasklist

    ! If not present PROBACTION_PERFORM_ALWAYS is assumed
    integer(i32), intent(in), optional :: iperformspec
!</input>
!</subroutine>

!!$    ! local variables
!!$    type(t_collection) :: rcollection
!!$    type(t_vectorBlock) :: rvectorTemp
!!$    type(t_vectorTask), pointer :: p_rtask
!!$    integer :: iperform
!!$
!!$    iperform = PROBACTION_PERFORM_ALWAYS
!!$    if (present(iperformSpec)) iperform = iperformSpec
!!$    
!!$    ! Iterate over all tasks
!!$    p_rtask => rtasklist
!!$    taskloop: do while (associated(p_rtask))
!!$
!!$      ! Do we have to perform this task?
!!$      if (iand(p_rtask%iperform, iperform) .eq. 0) then
!!$        p_rtask => p_rtask%p_rnextTask
!!$        cycle taskloop
!!$      end if
!!$      
!!$      !-------------------------------------------------------------------------
!!$      ! First part: create or duplicate vector
!!$      !-------------------------------------------------------------------------
!!$      select case(p_rtask%ctask)
!!$      case(PROBACTION_CREATE)
!!$
!!$        ! Do we have to assemble a scalar vector?
!!$        if (p_rtask%bisVectorScalar) then
!!$          
!!$          !---------------------------------------------------------------------
!!$          ! Create scalar vector
!!$          
!!$          ! Do we have spatial discretisations?
!!$          if (.not.associated(p_rtask%p_rspatialDiscr)) then
!!$            call output_line('Unable to find spatial discretisation',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$            call sys_halt()
!!$          end if
!!$          
!!$          ! Clear preexisting vector
!!$          call lsyssc_releaseVector(p_rtask%p_rvectorScalarDest)
!!$          
!!$          ! Create vector structure and set data type
!!$          call lsyssc_createVector(p_rtask%p_rspatialDiscr,&
!!$              p_rtask%p_rvectorScalarDest, p_rtask%nvar, .false., p_rtask%cdataType)
!!$
!!$        else
!!$
!!$          !---------------------------------------------------------------------
!!$          ! Create block vector
!!$
!!$          ! Do we have block discretisations?
!!$          if (.not.associated(p_rtask%p_rblockDiscr)) then
!!$            call output_line('Unable to find block discretisation',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$            call sys_halt()
!!$          end if
!!$
!!$          ! Clear preexisting block vector
!!$          call lsysbl_releaseVector(p_rtask%p_rvectorBlockDest)
!!$
!!$          ! Create block vector from discretisation
!!$          call lsysbl_createVecBlockByDiscr(p_rtask%p_rblockDiscr,&
!!$              p_rtask%p_rvectorBlockDest)
!!$          
!!$          ! Proceed with next task
!!$          p_rtask => p_rtask%p_rnextTask
!!$          
!!$          ! That`s it for block vectors
!!$          cycle taskloop
!!$        end if
!!$
!!$      case(PROBACTION_DUPLICATE)
!!$
!!$        !-----------------------------------------------------------------------
!!$        ! Duplicate scalar/block vector
!!$
!!$        if (p_rtask%bisVectorScalar) then
!!$          ! Duplicate scalar vector
!!$          call lsyssc_duplicateVector(p_rtask%p_rvectorScalarSrc,&
!!$              p_rtask%p_rvectorScalarDest, p_rtask%cdupStructure, p_rtask%cdupContent)
!!$        else
!!$          if (associated(p_rtask%p_rvectorScalarSrc)) then
!!$            ! Duplicate scalar vector as 1-block vector
!!$            call lsysbl_createVecFromScalar(p_rtask%p_rvectorScalarSrc, rvectorTemp)
!!$
!!$            ! Copy content of temporal vector by hand
!!$            p_rtask%p_rvectorBlockDest%NEQ         = rvectorTemp%NEQ
!!$            p_rtask%p_rvectorBlockDest%nblocks     = 1
!!$
!!$            allocate(p_rtask%p_rvectorBlockDest%RvectorBlock(1))
!!$            p_rtask%p_rvectorBlockDest%RvectorBlock(1) = rvectorTemp%RvectorBlock(1)
!!$
!!$            ! Release temporal vector
!!$            call lsysbl_releaseVector(rvectorTemp)
!!$
!!$          elseif(associated(p_rtask%p_rvectorBlockSrc)) then
!!$            ! Duplicate block vector
!!$            call lsysbl_duplicateVector(p_rtask%p_rvectorBlockSrc,&
!!$                p_rtask%p_rvectorBlockDest, p_rtask%cdupStructure, p_rtask%cdupContent)
!!$          else
!!$            call output_line('Neither scalar nor block vector can be duplicated',&
!!$                OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$            call sys_halt()
!!$          end if
!!$        end if
!!$
!!$      case default
!!$        call output_line('Unsupported action: '//sys_siL(p_rtask%ctask,3),&
!!$            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateMamtrix')
!!$        call sys_halt()
!!$      end select
!!$
!!$      !-------------------------------------------------------------------------
!!$      ! Second part: generate content (only for scalar vectors)
!!$      !-------------------------------------------------------------------------
!!$
!!$      ! What creation method is adopted?
!!$      select case (p_rtask%cmethod)
!!$      case(PROBMETH_VECTOR_EMPTY)
!!$
!!$        !---------------------------------------------------------------------
!!$        ! Create empty vector
!!$
!!$        call lsyssc_clearVector(p_rtask%p_rvectorScalarDest, 0.0_DP)
!!$
!!$      case(PROBMETH_VECTOR_UNITY)
!!$
!!$        !---------------------------------------------------------------------
!!$        ! Create unit vector
!!$
!!$        call lsyssc_clearVector(p_rtask%p_rvectorScalarDest, 1.0_DP)
!!$                
!!$      case(PROBMETH_VECTOR_LINF)
!!$
!!$        !---------------------------------------------------------------------
!!$        ! Create vector by bilinearform evaluation         
!!$
!!$        rcollection%p_rfparserQuickAccess1 => p_rtask%p_rfparser
!!$        if (associated(p_rtask%p_rcubatureInfo)) then
!!$          call linf_buildVectorScalar(p_rtask%p_rform, .true.,&
!!$              p_rtask%p_rvectorScalarDest, p_rtask%p_rcubatureInfo,&
!!$              problem_buildVectorSc_fparser_sim, rcollection)
!!$        elseif (associated(p_rtask%p_rspatialDiscr)) then
!!$          call linf_buildVectorScalar(p_rtask%p_rspatialDiscr,&
!!$              p_rtask%p_rform, .true., p_rtask%p_rvectorScalarDest,&
!!$              problem_buildVectorSc_fparser_sim, rcollection)
!!$        else
!!$          call output_line('Neither cubature info structure nor discretisation available',&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$          call sys_halt()
!!$        end if
!!$        
!!$      case(PROBMETH_VECTOR_VIRTUAL)
!!$        
!!$        !---------------------------------------------------------------------
!!$        ! Create virtual vector, aka, do nothing
!!$        
!!$      case default
!!$        call output_line('Unsupported method: '//sys_siL(p_rtask%cmethod,3),&
!!$            OU_CLASS_ERROR,OU_MODE_STD,'problem_updateVector')
!!$        call sys_halt()
!!$      end select
!!$      
!!$      ! Proceed with next task
!!$      p_rtask => p_rtask%p_rnextTask
!!$    end do taskloop
    
  end subroutine problem_updateAFCstab2

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_clearAFCstabTL(p_rtasklist)

!<description>
    ! This subroutine clears the stabilisation of AFC-type task list.
!</description>

!<input>
    ! task list
    type(t_afcstabTask), pointer :: p_rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_afcstabTask), pointer :: p_rtask,p_rtaskSave

    p_rtask => p_rtasklist
    do while(associated(p_rtask))
      p_rtaskSave => p_rtask
      p_rtask     => p_rtask%p_rnextTask
      deallocate(p_rtaskSave)
    end do
    
  end subroutine problem_clearAFCstabTL

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_infoAFCstabTL(rtasklist)

!<description>
    ! This subroutine outputs information about the task list of
    ! stabilisation structures of AFC-type.
!</description>

!<input>
    ! task list
    type(t_afcstabTask), intent(in), target :: rtasklist
!</input>
!</subroutine>

    ! local variables
    type(t_afcstabTask), pointer :: p_rtask

    p_rtask => rtasklist
    do while(associated(p_rtask))

      select case(p_rtask%ctask)
      case (PROBACTION_NONE)
        call output_lbrk()
        call output_line ('AFCstabTask: DO NOTHING')

      case (PROBACTION_CREATE)
        call output_lbrk()
        call output_line ('AFCstabTask: CREATE')

      case (PROBACTION_DUPLICATE)
        call output_lbrk()
        call output_line ('AFCstabTask: DUPLICATE')

      case default
        call output_line ('AFCstabTask: UNSUPPORTED')
      end select

!!$      call output_line ('-----------------------')
!!$      call output_line ('vector type               : '//merge('SCALAR', 'BLOCK ', p_rtask%bisVectorScalar))
!!$      call output_line ('section name              : '//trim(p_rtask%ssectionName))
!!$      call output_line ('destination problem       : '//trim(p_rtask%p_rproblemDest%cproblem))
!!$      call output_line ('destination problem level : '//trim(sys_siL(p_rtask%p_rproblemLevelDest%ilev,15)))
!!$      call output_line ('destination vector scalar : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorScalarDest)))
!!$      call output_line ('destination vector block  : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorBlockDest)))
!!$      if (associated(p_rtask%p_rproblemSrc))&
!!$      call output_line ('source problem            : '//trim(p_rtask%p_rproblemSrc%cproblem))
!!$      if (associated(p_rtask%p_rproblemLevelSrc))&
!!$      call output_line ('source problem level      : '//trim(sys_siL(p_rtask%p_rproblemLevelSrc%ilev,15)))
!!$      call output_line ('source vector scalar      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorScalarSrc)))
!!$      call output_line ('source vector block       : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rvectorBlockSrc)))
!!$      call output_line ('iperform                  : '//trim(sys_siL(int(p_rtask%iperform),15)))
!!$      call output_line ('method                    : '//trim(sys_siL(p_rtask%cmethod,15)))
!!$      call output_line ('data type                 : '//trim(sys_siL(p_rtask%cdatatype,15)))
!!$      call output_line ('nvar                      : '//trim(sys_siL(p_rtask%nvar,15)))
!!$      call output_line ('spatial discretisation:     '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rspatialDiscr)))
!!$      call output_line ('block discretisation      : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rblockDiscr)))
!!$      call output_line ('cubature info structure   : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rcubatureInfo)))
!!$      call output_line ('function parser           : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rfparser)))
!!$      call output_line ('linear form evaluator     : '//merge('ASSOCIATED    ','NOT ASSOCIATED',&
!!$                                                              associated(p_rtask%p_rform)))

      p_rtask => p_rtask%p_rnextTask
    end do
    
  end subroutine problem_infoAFCstabTL

  

end module problem
