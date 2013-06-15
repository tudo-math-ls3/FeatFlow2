!##############################################################################
!# ****************************************************************************
!# <name> collection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a dynamic list of named values in the memory.
!# Values can be added to the list and changed if necessary.
!# This corresponds to the assignment "NAME = VALUE" with VALUE being
!# - character
!# - Character string
!# - double real
!# - integer
!# - a vector - scalar and block type
!# - a matrix - scalar and block type
!# - triangulation structures
!# - discrete BC structures
!# - geometry information structures
!# - linear solver configurations
!# - a parameter list
!# - a scalar discretisation structure
!# - a block discretisation structure
!# - a parser structure
!# - another collection
!# - small, fixed-size integer arrays
!# - small, fixed-size double precision arrays
!# - a grid adaptation structure
!# - a stabilisation structure
!# - a linked list
!# - an arraylist
!# - a map
!# - a graph structure
!# - a group finite element set and block type
!#
!# The list supports 'level tags' and 'section tags' to group values:
!#
!# 1.) There is one unnamed section ('') in the collection for general purpose
!#     variables. All other sections are named and independent from each other.
!#     Sections can be added if necessary.
!# 2.) Each section contains of a 'level 0' list for general, level-independent
!#     data. Furthermore it contains a level structure for level 1..n
!#     where n can by increased dynamically. All level are independent from
!#     each other.
!#
!# Names of sections and variables are case-insensitive. Allowed characters
!# for names are: 'A'..'Z', '0'..'9', '-', '_', '(', ')'.
!#
!# A collection is similar to a large INI file containing all basic
!# information. One can think of the following structure:
!#
!# <verb>
!# -------------------SNIP--------------------
!#
!# PROBLEM_TYPE   = 'CC2DSTATIONARY'
!# DIMENSION      = 2
!# NLMIN          = 1
!# NLMAX          = 3
!#
!# [CC2DSTATIONARY]
!# SOLVER         = 'NONLINEARDEFCOR'
!# GRIDADAPT      = 'GRIDADAPT'
!# ...
!# GEOMETRY       = (geometry structure)
!# LINSOL         = (linear solver structure for all levels)
!#
!#   [LEVEL1-Substructure]
!#   TRIANGULATION  = (triangulation structure on level 1)
!#   DISCRETISATION = (discr. structure on level 1)
!#   DISCRETEBC     = (discrete BC structure on level 1)
!#   MASSMATRIX     = (matrix structure for mass matrix on level 1)
!#   LAPLACEMATRIX  = (matrix structrue for laplace matrix on level 1)
!#   SYSTEMMATRIX   = (matrix structrue for system matrix on level 1)
!#   SOLUTION       = (vector structure for solution vector on level 1)
!#   RHS            = (vector structure for RHS vector on level 1)
!#   ...
!#
!#   [LEVEL2-Substructure]
!#   TRIANGULATION  = (triangulation structure on level 2)
!#   DISCRETISATION = (discr. structure on level 2)
!#   DISCRETEBC     = (discrete BC structure on level 2)
!#   MASSMATRIX     = (matrix structure for mass matrix on level 2)
!#   LAPLACEMATRIX  = (matrix structrue for laplace matrix on level 2)
!#   SYSTEMMATRIX   = (matrix structrue for system matrix on level 2)
!#   SOLUTION       = (vector structure for solution vector on level 2)
!#   RHS            = (vector structure for RHS vector on level 2)
!#   ...
!#
!#   [LEVEL3-Substructure]
!#   TRIANGULATION  = (triangulation structure on level 3)
!#   DISCRETISATION = (discr. structure on level 3)
!#   DISCRETEBC     = (discrete BC structure on level 3)
!#   MASSMATRIX     = (matrix structure for mass matrix on level 3)
!#   LAPLACEMATRIX  = (matrix structrue for laplace matrix on level 3)
!#   SYSTEMMATRIX   = (matrix structrue for system matrix on level 3)
!#   SOLUTION       = (vector structure for solution vector on level 3)
!#   RHS            = (vector structure for RHS vector on level 3)
!#   ...
!#
!# [NONLINEARDEFCOR]
!# ... (parameters for the nonlinear defect correction)
!#
!# [GRIDADAPT]
!# ... (parameters for the grid adaption)
!#
!# -------------------SNIP--------------------
!# </verb>
!# Additionally to these 'named variables', the collection contains
!# a small integer and double precision 'temporary quick access'
!# array directly in the t_collection structure with size of at
!# least 128 elemens.
!# Care must be taken when an application wants to use these arays.
!# Not maintained by the collection internally, they are of
!# *pure temporary nature*. Long-life information should never be stored
!# in there! They can be used e.g. in assembly routines to provide
!# callback routines with information that must quickly be accessed
!# without the need of asking the collection for a named variable.
!# Their content can be described at best as `only valid for one
!# operation` (e.g. for one matrix asembly).
!#
!# The following routines can be used to maintain a parameter
!# list:
!# 1.) collct_init
!#     -> Initialises an empty list of named scalar values.
!#
!# 2.) collct_done
!#     -> Cleans up a list, releases all allocated memory from the heap.
!#
!# 3.) collct_getmaxlevel
!#     -> Determines the number of levels.
!#
!# 4.) collct_addlevel
!#     -> Adds a new level.
!#
!# 5.) collct_deletelevel
!#     -> Deletes the maximum level
!#        (not yet implemented)
!#
!# 6.) collct_queryvalue
!#     -> Determines whether or not a parameter exists
!#
!# 7.) collct_gettype
!#     -> Get the type of a variable.
!#
!# 8.) collct_getvalue_xxx
!#     -> Get a value from the list. xxx corresponds to the type of
!#        the variable.
!#
!# 9.) collct_setvalue_xxx
!#     -> Modifies a value in the list. Add a value if it does not exist.
!#        xxx corresponds to the type of  the variable.
!#
!# 10.) collct_deletevalue
!#      -> Deletes a value from the list.
!#
!# 11.) collct_printStatistics
!#      -> Prints out the current content of the structure to the terminal
!#
!# 12.) collct_addsection
!#      -> Adds a section to the collection
!#
!# 13.) collct_deletesection
!#      -> Removes a section from the collection, deletes all its content
!#
!# 14.) collct_settag
!#      -> Sets the user defined tag of a value
!#
!# 15.) collct_gettag
!#      -> Returns the user defined tag of a value
!#
!# 16.) collct_copyQuickAccess
!#      -> Copies all quick access arrays from one collection structure
!#         to another collection structure
!#
!# Auxiliary routines:
!#
!#  1.) collct_cleanupvalue
!#      -> Cleans up a value structure.
!#
!# GUILDLINE \\
!# --------- \\
!# Here a small guidline how and when to use a collection:
!#
!# HOLDING PROBLEM DEPENDENT DATA FOR CALLBACK ROUTINES,
!# PRINCIPLE OF MINIMALISM
!#
!#  Use a collection to store data that must be available to callback routines
!#  for solving a problem. Data that is only used in the main application
!#  should be stored in a problem-dependent structure.
!#  Try to save as few information in it as possible. This helps you not to
!#  loose control about what is in the container and what not.
!#
!# NO PROBLEM-SPECIFIC STRUCTURES
!#
!#  The collection is a black-box container. It is designed to hold low- and
!#  mid-level information (numbers, strings, as well as linear algebra stuff
!#  like vectors, matrices) that is used for solving linear problems. It is
!#  clearly NOT designed that a user adds problem-specific data structures
!#  there. Trying this will result in circular dependencies! Adding information
!#  from a problem-dependent structure to the collection with the 'add'
!#  routines is completely right, but do not extent the collection
!#  by new data types / structures unless you are absolutely sure that this
!#  is relevant for all users and that you will not get circular dependencies!
!#
!# USE THE COLLECTION AS AN INDEX
!#
!#  Although it is possible to use the collection as a data container for most
!#  (if not all) of your information, you should not use the collection
!#  as the main data container. Put your problem-dependent data into
!#  a problem-dependent structure and add pointers to them to the collection
!#  (think of pointers to vectors/matrices)
!#  Therefore, use the collection as an *index* to your data and not as
!#  the main data container. This helps you to prevent memory holes,
!#  since the collection does not delete anything from the heap when you
!#  remove a variable from it - it simply deletes a pointer!
!#
!# USE THE TEMPORARY QUICK ACCESS ARRAYS FOR INFORMATION TO BE ACCESSED QUICKLY
!#
!#  Although the quick-access arrays are only of temporary nature,
!#  you can use them to store information of simple type that you often
!#  need in a callback routine. Make use of this feature to get speed,
!#  but do not use this feature for storing long-life information!
!#  Example:
!#
!# <code>
!#    PROGRAM Test
!#      TYPE(t_collection) :: rcollection
!#      ...
!#      rcollection%Iquickaccess(1) = 1         ! Number of boundary component e.g.
!#      rcollection%Iquickaccess(2) = 3         ! Number of Dirichlet-segment e.g.
!#      rcollection%Dquickaccess(1) = 1.0_DP    ! Inflow speed e.g.
!#      rcollection%SquickAccess(1) = 'PARPROF' ! Name of the profile
!#      CALL bcasm_newDirichletBConRealBD (...,callback,rcollection)
!#      ...
!#    END PROGRAM
!#
!#    SUBROUTINE callback(...,rbcRegion,...,rcollection,...,Dvalues)
!#      ...
!#      iboundaryComponent = rcollection%Iquickaccess(1)
!#      iinflowSegment = rcollection%Iquickaccess(2)
!#      dinfowSpeed    = rcollection%Dquickaccess(1)
!#      sname          = rcollection%Squickaccess(1)
!#      IF ((iboundaryComponent .EQ. rbcRegion%iboundCompIdx) .AND.
!#          (iinflowSegment .EQ. rbcRegion%iboundSegIdx)) THEN
!#        Dvalues(1) = dinflowSpeed
!#      END IF
!#      ...
!#    END SUBROUTINE
!# </code>
!#
!#  rcollection%p_rvectorQuickAccessX and rcollection%p_rmatrixQuickAccessX
!#  provides pointers for some user defined matrices and vectors which
!#  can be accessed directly.
!#
!#  The rcollection%p_rnextCollection pointer is also user defined. Here,
!#  the user can set a pointer to another collection structure if necessary.
!#  This allows to build a small list of collection structures that can
!#  directly be accessed.
!#
!# STORE ONLY SMALL-SIZE ARRAYS
!#
!#  It is possible to store integer / real arrays in the collection. These
!#  should only be small and the application must keep track of the size
!#  of these arrays. To store larger arrays, use the STORAGE module to
!#  allocate memory and store the handle of the memory block in the collection!
!#
!# </purpose>
!##############################################################################

module collection

!$use omp_lib
  use fsystem
  use io
  use genoutput
  use paramlist, only: t_parlist
  use boundary, only: t_boundary, t_boundaryRegion
  use boundarycondition, only: t_boundaryConditions
  use triangulation, only: t_triangulation
  use spatialdiscretisation, only: t_spatialDiscretisation, &
                                   t_blockDiscretisation,t_scalarCubatureInfo
  use discretebc, only: t_discreteBC
  use linearsystemscalar, only: t_vectorScalar, t_matrixScalar
  use linearsystemblock, only: t_vectorBlock, t_matrixBlock
  use linearsolver, only: t_linsolNode
  use multilevelprojection, only: t_interlevelProjectionScalar, &
                                  t_interlevelProjectionBlock
  use filtersupport, only: t_filterChain
  use fparser, only: t_fparser
  use geometry, only: t_geometryObject, t_particleCollection, t_particleCollection3D
  use hadaptaux, only: t_hadapt
  use afcstabbase, only: t_afcstab
  use statistics, only: t_timer
  use graph, only: t_graph
  use meshhierarchy, only: t_meshHierarchy
  use timescalehierarchy, only: t_timescaleHierarchy
  use fespacehierarchybase, only: t_feHierarchy, t_feSpaceLevel
  use multilevelprojection, only: t_interlevelProjectionHier
  use groupfembase, only: t_groupFEMSet, t_groupFEMBlock

  implicit none

  private

!<constants>

!<constantblock>

  ! Maximum length of a section name.
  integer, parameter, public :: COLLCT_MLSECTION = 64

  ! Maximum length of the name of a value: 32 characters
  integer, parameter, public :: COLLCT_MLNAME = 32

  ! Minimum number of free 'slots' per 'level' subgroup.
  ! If there are too many values in a 'level' subgroup, the
  ! structure is dynamically extended in terms of collct_NVALUES
  ! entries.
  integer, parameter, public :: COLLCT_NVALUES = 32

  ! Minimum number of level subgroups (additionally to the 0th).
  ! If there are too many levels in a parameter block, the
  ! structure is dynamically extended in terms of COLLCT_NLEVELS
  ! entries.
  integer, parameter, public :: COLLCT_NLEVELS = 9

  ! Minimum number of sections - additionally to the unnamed one.
  ! If there are too many sections in a collection block, the
  ! structure is dynamically extended in terms of COLLCT_NSECTIONS
  ! entries.
  integer, parameter, public :: COLLCT_NSECTIONS = 5

  ! Length of each 'quick-access' array in the collection.
  integer, parameter, public :: COLLCT_QALENGTH = 128

  ! Number of 'quick-access' strings in the collection
  integer, parameter, public :: COLLCT_QASTRINGS = 4
!</constantblock>

!<constantblock description="Type identifier for values">

  ! Undefined value
  integer, parameter, public :: COLLCT_UNDEFINED    = 0

  ! Character value
  integer, parameter, public :: COLLCT_CHARACTER    = 1

  ! String value
  integer, parameter, public :: COLLCT_STRING       = 2

  ! Integer value
  integer, parameter, public :: COLLCT_INTEGER      = 3

  ! Double precision REAL value
  integer, parameter, public :: COLLCT_REAL         = 4

  ! Scalar discretisation structure
  integer, parameter, public :: COLLCT_DISCR        = 5

  ! Block discretisation structure
  integer, parameter, public :: COLLCT_BLDISCR      = 5

  ! Triangulation structure
  integer, parameter, public :: COLLCT_TRIA         = 7

  ! A scalar vector
  integer, parameter, public :: COLLCT_SCAVECTOR    = 8

  ! A scalar matrix
  integer, parameter, public :: COLLCT_SCAMATRIX    = 9

  ! A block vector
  integer, parameter, public :: COLLCT_BLKVECTOR    = 10

  ! A block matrix
  integer, parameter, public :: COLLCT_BLKMATRIX    = 11

  ! A parameter structure
  integer, parameter, public :: COLLCT_PARAMETERS   = 12

  ! A linear solver structure
  integer, parameter, public :: COLLCT_LINSOL       = 13

  ! Discretised-boundary-conditions structure
  integer, parameter, public :: COLLCT_DISCRBC      = 14

  ! A domain
  integer, parameter, public :: COLLCT_BOUNDARY     = 15

  ! Scalar analytic boundary conditions
  integer, parameter, public :: COLLCT_BOUNDARYCOND = 16

  ! Scalar interlevel projection structure
  integer, parameter, public :: COLLCT_INTERLVPRJSC = 17

  ! Block interlevel projection structure
  integer, parameter, public :: COLLCT_INTERLVPRJ   = 18

  ! A parser structure for a 'precompiled' expression
  integer, parameter, public :: COLLCT_FPARSER      = 19

  ! A filter chain
  integer, parameter, public :: COLLCT_FILTERCHAIN  = 20

  ! Integer value type array
  integer, parameter, public :: COLLCT_INTEGERARR   = 21

  ! Double precision REAL value type array
  integer, parameter, public :: COLLCT_REALARR      = 22

  ! Geometry object
  integer, parameter, public :: COLLCT_GEOMETRY     = 24

  ! The collection structure itself
  integer, parameter, public :: COLLCT_COLLECTION   = 25

  ! Grid adaptation object
  integer, parameter, public :: COLLCT_HADAPT       = 26

  ! Stabilisation structure
  integer, parameter, public :: COLLCT_AFCSTAB      = 27

  ! Timer structure
  integer, parameter, public :: COLLCT_TIMER        = 28

  ! Graph structure
  integer, parameter, public :: COLLCT_GRAPH        = 29

  ! Boundary region
  integer, parameter, public :: COLLCT_BDREGION     = 30

  ! Particles structure
  integer, parameter, public :: COLLCT_PARTICLES    = 31

  ! Particles structure3D
  integer, parameter, public :: COLLCT_PARTICLES3D    = 32

  ! Mesh hierarchy
  integer, parameter, public :: COLLCT_MSHHIERARCHY   = 33

  ! FE space
  integer, parameter, public :: COLLCT_FESPACE        = 34

  ! FE space hierarchy
  integer, parameter, public :: COLLCT_FEHIERARCHY    = 35

  ! Time scale hierarchy
  integer, parameter, public :: COLLCT_TSHIERARCHY    = 36

  ! Multilevel projection hierarchy
  integer, parameter, public :: COLLCT_MLPRJHIERARCHY = 37

  ! Group finite element set
  integer, parameter, public :: COLLCT_GROUPFEMSET    = 38

  ! Group finite element block
  integer, parameter, public :: COLLCT_GROUPFEMBLOCK  = 39

  ! Cubature information structure
  integer, parameter, public :: COLLCT_CUBINFO        = 40

  ! Generic object
  integer, parameter, public :: COLLCT_GENOBJ          = 41

  ! Pointer to Integer value type array
  integer, parameter, public :: COLLCT_INTEGERARRP     = 42

  ! Pointer to Double precision REAL value type array
  integer, parameter, public :: COLLCT_REALARRP        = 43

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This structure realises a value in the scalar collection.
  ! It contains a type and some data variables that contain the value.

  type t_collctValue

    private

    ! The type of the value (COLLCT_UNDEFINED, COLLCT_CHARACTER,...)
    integer :: itype = COLLCT_UNDEFINED

    ! A user defined tag. Standard value is =0.
    integer :: itag = 0

    ! Name of the value
    character(LEN=COLLCT_MLNAME) :: sname

    ! Character value
    character :: csvalue = " "

    ! String value
    character, dimension(:), pointer :: p_svalue => null()

    ! Integer value
    integer :: ivalue = 0

    ! Double precision value
    real(DP) :: dvalue = 0.0_DP

    ! Pointer to a spatial discretisation structure
    type(t_spatialDiscretisation), pointer      :: p_rdiscretisation => null()

    ! Pointer to a block discretisation structure
    type(t_blockDiscretisation), pointer        :: p_rblockDiscretisation => null()

    ! Pointer to a triangulation structure
    type(t_triangulation), pointer              :: p_rtriangulation => null()

    ! Pointer to a scalar vector
    type(t_vectorScalar), pointer               :: p_rvectorScalar => null()

    ! Pointer to a scalar matrix
    type(t_matrixScalar), pointer               :: p_rmatrixScalar => null()

    ! Pointer to a block vector
    type(t_vectorBlock), pointer                :: p_rvector => null()

    ! Pointer to a block matrix
    type(t_matrixBlock), pointer                :: p_rmatrix => null()

    ! Pointer to a parameter list
    type(t_parlist), pointer                    :: p_rparlist => null()

    ! Pointer to a linear solver structure
    type(t_linsolNode), pointer                 :: p_rlinearSolver => null()

    ! A structure containing discretised boundary conditions
    type(t_discreteBC), pointer                 :: p_rdiscreteBC     => null()

    ! Pointer to a domain
    type(t_boundary), pointer                   :: p_rboundary => null()

    ! Pointer to a boundary region
    type(t_boundaryRegion), pointer             :: p_rboundaryRegion => null()

    ! Pointer to a geometry object
    type(t_geometryObject), pointer             :: p_rgeometry => null()

    ! Pointer to a particle collection
    type(t_particleCollection), pointer         :: p_rparticles => null()

    ! Pointer to a particle collection
    type(t_particleCollection3D), pointer       :: p_rparticles3D => null()

    ! Pointer to scalar boundary conditions
    type(t_boundaryConditions), pointer         :: p_rboundaryConditions => null()

    ! Pointer to a scalar interlevel projection structure
    type(t_interlevelProjectionScalar), pointer :: p_rilvprojectionSc => null()

    ! Pointer to a scalar interlevel projection structure
    type(t_interlevelProjectionBlock), pointer  :: p_rilvprojection => null()

    ! Pointer to a parser object
    type(t_fparser), pointer                    :: p_rparser => null()

    ! Pointer to a filter chain
    type(t_filterChain), dimension(:), pointer  :: p_RfilterChain => null()

    ! Pointer to an integer precision array
    integer, dimension(:), pointer              :: p_Iarray => null()

    ! Pointer to an real precision array
    real(DP), dimension(:), pointer             :: p_Darray => null()

    ! Pointer to a collection structure
    type(t_collection), pointer                 :: p_rcollection => null()

    ! Pointer to a grid adaptation structure
    type(t_hadapt), pointer                     :: p_rhadapt => null()

    ! Pointer to a stabilisation structure
    type(t_afcstab), pointer                    :: p_rafcstab => null()

    ! Pointer to a timer structure
    type(t_timer), pointer                      :: p_rtimer => null()

    ! Pointer to a graph structure
    type(t_graph), pointer                      :: p_rgraph => null()

    ! Pointer to a mesh hierarchy structure
    type(t_meshHierarchy), pointer              :: p_rmeshHierarchy => null()

    ! Pointer to a FE space hierarchy structure
    type(t_feHierarchy), pointer                :: p_rfeHierarchy => null()

    ! Pointer to a FE space hierarchy structure
    type(t_feSpaceLevel), pointer               :: p_rfeSpace => null()

    ! Pointer to a time scale hierarchy structure
    type(t_timescaleHierarchy), pointer         :: p_rtsHierarchy => null()

    ! Pointer to a multilevel projection hierarchy structure
    type(t_interlevelProjectionHier), pointer   :: p_rmlprjHierarchy => null()

    ! Pointer to a group finite element set
    type(t_groupFEMSet), pointer                :: p_rgroupFEMSet => null()

    ! Pointer to a group finite element block
    type(t_groupFEMBlock), pointer              :: p_rgroupFEMBlock => null()

    ! Pointer to a cubature information structure
    type(t_scalarCubatureInfo), pointer         :: p_rcubatureInfo => null()

    ! Pointer to a generic object structure
    type(t_genericObject), pointer              :: p_rgenericObject => null()
  end type

  public :: t_collctValue

!</typeblock>

!<typeblock>

  ! This structure realises a level structure containing multiple values.

  type t_collctLevel

    private

    ! Actual number of values in this section.
    integer :: ivalueCount = 0

    ! Indicates whether the first ivalueCount values are 'en-block'
    ! or if there are 'holes' inbetween because of deleted values.
    logical :: bisFull = .true.

    ! A list of values.
    type(t_collctValue), dimension(:), pointer :: p_Rvalues => null()

  end type

  public :: t_collctLevel

!</typeblock>

!<typeblock>

  ! This structure realises a section. Each section contains a list of
  ! level substructures and one substructure not assigned to a level,
  ! containing level independent parameters.

  type t_collctSection

    private

    ! The name of the section. '' identifies an/the unnamed section.
    character(LEN=COLLCT_MLSECTION) :: ssectionName = ''

    ! Actual number of levels in this section. There is at least
    ! one section - the unnamed section. If this value is =0, the parameter
    ! list is currently not initialised.
    integer :: ilevelCount = 0

    ! The 0th level collecting general variables, level-independent
    type(t_collctLevel) :: rlevel0

    ! A list of levels or NULL(), if there are no additional 'level'
    ! subgroups than level 0.
    type(t_collctLevel), dimension(:), pointer :: p_Rlevels => null()

  end type

  public :: t_collctSection

!</typeblock>

!<typeblock>

  ! This structure realises a collection.
  ! A collection contains one unnamed section (containing general data)
  ! and arbitrary many named sections.
  ! Each section contains a level structure for level 0 (containing
  ! level independent data) and a number of levels where data can
  ! be stored.

  type t_collection

    ! USER DEFINED:
    ! Quick-access array for integers. This is a short, temporary
    ! array that can be used by the application to save intermediate
    ! values (e.g. some information that is to pass and to be accessed
    ! by callback routines).
    ! No long life information should be stored here. Long-life
    ! information should be named and added as such to the collection
    ! directly. The collection itself does not maintain this array!
    integer, dimension(COLLCT_QALENGTH) :: IquickAccess = 0

    ! USER DEFINED:
    ! Quick-access array for reals. This is a short, temporary
    ! array that can be used by the application to save intermediate
    ! values (e.g. some information that is to pass and to be accessed
    ! by callback routines).
    ! No long life information should be stored here. Long-life
    ! information should be named and added as such to the collection
    ! directly. The collection itself does not maintain this array!
    real(DP), dimension(COLLCT_QALENGTH) :: DquickAccess = 0.0_DP

    ! USER DEFINED:
    ! A simple string array for direct access.
    character(LEN=SYS_STRLEN), dimension(COLLCT_QASTRINGS) :: SquickAccess

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess1 => null()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess2 => null()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess3 => null()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess4 => null()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess5 => null()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess6 => null()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess7 => null()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    type(t_vectorBlock), pointer :: p_rvectorQuickAccess8 => null()

    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    type(t_matrixBlock), pointer :: p_rmatrixQuickAccess1 => null()

    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    type(t_matrixBlock), pointer :: p_rmatrixQuickAccess2 => null()

    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    type(t_matrixBlock), pointer :: p_rmatrixQuickAccess3 => null()

    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    type(t_matrixBlock), pointer :: p_rmatrixQuickAccess4 => null()

    ! USER DEFINED:
    ! A quick access pointer to a parameter list.
    type(t_parlist), pointer :: p_rparlistQuickAccess1 => null()

    ! USER DEFINED:
    ! A quick access pointer to a parameter list.
    type(t_parlist), pointer :: p_rparlistQuickAccess2 => null()

    ! USER DEFINED:
    ! A quick access pointer to a parameter list.
    type(t_parlist), pointer :: p_rparlistQuickAccess3 => null()

    ! USER DEFINED:
    ! A quick access pointer to a parameter list.
    type(t_parlist), pointer :: p_rparlistQuickAccess4 => null()

    ! USER DEFINED:
    ! A quick acces pointer to a function parser
    type(t_fparser), pointer :: p_rfparserQuickAccess1 => null()

    ! USER DEFINED:
    ! A quick acces pointer to a function parser
    type(t_fparser), pointer :: p_rfparserQuickAccess2 => null()

    ! USER DEFINED:
    ! A quick acces pointer to a function parser
    type(t_fparser), pointer :: p_rfparserQuickAccess3 => null()

    ! USER DEFINED:
    ! A quick acces pointer to a function parser
    type(t_fparser), pointer :: p_rfparserQuickAccess4 => null()

    ! USER DEFINED:
    ! A quick access pointer to another collection structure.
    ! Can be used to link collection structures.
    type(t_collection), pointer :: p_rnextCollection => null()

    ! Actual number of sections in this collection.
    ! This is at least one, as every collection contains at least an
    ! unnamed section.
    integer :: isectionCount = 0

    ! A list of levels or NULL(), if there are no additional 'level'
    ! subgroups than level 0.
    type(t_collctSection), dimension(:), pointer :: p_Rsections => null()

  end type

  public :: t_collection

!</typeblock>

!</types>

  interface collct_addlevel
    module procedure collct_addlevel_indir
    module procedure collct_addlevel_direct
    module procedure collct_addlevel_all
  end interface

  interface collct_queryvalue
    module procedure collct_queryvalue_direct
    module procedure collct_queryvalue_indir
  end interface

  interface collct_getmaxlevel
    module procedure collct_getmaxlevel_direct
    module procedure collct_getmaxlevel_indir
  end interface

  public :: collct_init
  public :: collct_done
  public :: collct_addlevel
  public :: collct_addsection
  public :: collct_deletesection
  public :: collct_cleanupvalue
  public :: collct_deletevalue
  public :: collct_printStatistics
  public :: collct_copyQuickAccess

  public :: collct_setvalue_char
  public :: collct_setvalue_string
  public :: collct_setvalue_int
  public :: collct_setvalue_intarr
  public :: collct_setvalue_iarrp
  public :: collct_setvalue_real
  public :: collct_setvalue_realarr
  public :: collct_setvalue_rarrp
  public :: collct_setvalue_discr
  public :: collct_setvalue_bldiscr
  public :: collct_setvalue_tria
  public :: collct_setvalue_bdry
  public :: collct_setvalue_bdreg
  public :: collct_setvalue_bc
  public :: collct_setvalue_vecsca
  public :: collct_setvalue_matsca
  public :: collct_setvalue_vec
  public :: collct_setvalue_mat
  public :: collct_setvalue_parlst
  public :: collct_setvalue_linsol
  public :: collct_setvalue_discbc
  public :: collct_setvalue_ilvpsc
  public :: collct_setvalue_ilvp
  public :: collct_setvalue_coll
  public :: collct_setvalue_pars
  public :: collct_setvalue_geom
  public :: collct_setvalue_fchn
  public :: collct_setvalue_hadapt
  public :: collct_setvalue_afcstab
  public :: collct_setvalue_timer
  public :: collct_setvalue_mshh
  public :: collct_setvalue_fesp
  public :: collct_setvalue_feh
  public :: collct_setvalue_tsh
  public :: collct_setvalue_mlprjh
  public :: collct_setvalue_graph
  public :: collct_setvalue_particles
  public :: collct_setvalue_particles3D
  public :: collct_setvalue_gfemset
  public :: collct_setvalue_gfemblk
  public :: collct_setvalue_cubinfo
  public :: collct_setvalue_genobj

  public :: collct_getmaxlevel
  public :: collct_getmaxlevel_direct
  public :: collct_getmaxlevel_indir
  public :: collct_queryvalue
  public :: collct_queryvalue_indir
  public :: collct_queryvalue_direct
  public :: collct_gettype
  public :: collct_settag
  public :: collct_gettag

  public :: collct_getvalue_struc
  public :: collct_getvalue_string
  public :: collct_getvalue_intarr
  public :: collct_getvalue_realarr
  public :: collct_getvalue_iarrp
  public :: collct_getvalue_rarrp
  public :: collct_getvalue_char
  public :: collct_getvalue_int
  public :: collct_getvalue_real
  public :: collct_getvalue_discr
  public :: collct_getvalue_bldiscr
  public :: collct_getvalue_tria
  public :: collct_getvalue_bdry
  public :: collct_getvalue_bdreg
  public :: collct_getvalue_bc
  public :: collct_getvalue_vecsca
  public :: collct_getvalue_matsca
  public :: collct_getvalue_vec
  public :: collct_getvalue_mat
  public :: collct_getvalue_linsol
  public :: collct_getvalue_discbc
  public :: collct_getvalue_parlst
  public :: collct_getvalue_ilvpsc
  public :: collct_getvalue_ilvp
  public :: collct_getvalue_coll
  public :: collct_getvalue_pars
  public :: collct_getvalue_geom
  public :: collct_getvalue_particles
  public :: collct_getvalue_particles3D
  public :: collct_getvalue_fchn
  public :: collct_getvalue_hadapt
  public :: collct_getvalue_afcstab
  public :: collct_getvalue_timer
  public :: collct_getvalue_graph
  public :: collct_getvalue_mshh
  public :: collct_getvalue_fesp
  public :: collct_getvalue_feh
  public :: collct_getvalue_tsh
  public :: collct_getvalue_mlprjh
  public :: collct_getvalue_gfemset
  public :: collct_getvalue_gfemblk
  public :: collct_getvalue_cubinfo
  public :: collct_getvalue_genobj

contains

  ! ***************************************************************************

!<subroutine>

  subroutine collct_initlevel (rcollctLevel)

!<description>

    ! Internal subroutine: Initialise a newly created level subgroup.

!</description>

!<inputoutput>

    ! The section to initialise
    type(t_collctLevel), intent(inout) :: rcollctLevel

!</inputoutput>
!</subroutine>

    ! Nullify the level pointer. Is allocated when the first entry is put
    ! into the list.
    nullify(rcollctLevel%p_Rvalues)

    ! No variables in here for the moment
    rcollctLevel%ivalueCount = 0

    ! No holes.
    rcollctLevel%bisFull = .true.

  end subroutine collct_initlevel

  ! ***************************************************************************

!<subroutine>

  subroutine collct_donelevel (rcollctLevel)

!<description>

  ! Internal subroutine: Releases a level and all values inside from memory.

!</description>

!<inputoutput>

    ! The section to release
    type(t_collctLevel), intent(inout) :: rcollctLevel

!</inputoutput>
!</subroutine>

    ! local variable
    integer :: i

    ! Deallocate the value list on the current level if there is one.
    if (associated(rcollctLevel%p_Rvalues)) then
      ! Release and deallocate
      do i = 1, size(rcollctLevel%p_Rvalues)
        call collct_cleanupvalue(rcollctLevel%p_Rvalues(i))
      end do
      deallocate(rcollctLevel%p_Rvalues)
    end if

    ! No variables in here for the moment
    rcollctLevel%ivalueCount = 0

    ! No holes anymore.
    rcollctLevel%bisFull = .true.

  end subroutine collct_donelevel

  ! ***************************************************************************

!<subroutine>

  subroutine collct_realloclevel (rcollctLevel, inewsize)

!<description>

    ! Internal subroutine: Reallocate a level.
    ! This increases the size of a level subgroup by
    ! reallocation of the arrays.

!</description>

!<input>

    ! The new 'size' of the section, i.e. the new number of parameters,
    ! the level should be able to handle.
    integer, intent(in) :: inewsize

!</input>

!<inputoutput>

    ! The section to reallocate.
    type(t_collctLevel), intent(inout) :: rcollctLevel

!</inputoutput>
!</subroutine>

    ! Pointers to new lists for replacing the old.
    type(t_collctValue), dimension(:), pointer :: p_Rvalues

    ! local variables
    integer :: sz

    sz = max(size(rcollctLevel%p_Rvalues), inewsize)

    if (size(rcollctLevel%p_Rvalues) .eq. sz) return ! nothing to do

    ! Allocate the pointers for the new list
    allocate(p_Rvalues(sz))

    ! Copy the content of the old ones
    p_Rvalues(1:sz) = rcollctLevel%p_Rvalues (1:sz)

    ! Throw away the old array, replace by the new one
    deallocate(rcollctLevel%p_Rvalues)

    rcollctLevel%p_Rvalues => p_Rvalues

  end subroutine collct_realloclevel

  ! ***************************************************************************

!<subroutine>

  subroutine collct_addvalue (rcollection, ssectionName, sparamName, itype, &
                              ilevel, p_rvalue)

!<description>

    ! Internal subroutine: Add a new value structure to a level and return
    ! a pointer to it. Returns NULL() if the section or the level does not
    ! exist.

!</description>

!<input>

    ! The name of the section; '' identifies the unnamed section
    character(LEN=*), intent(in) :: ssectionName

    ! The name of the value to add. Must be <> ''!
    character(LEN=*), intent(in) :: sparamName

    ! The type of the parameter
    integer, intent(in) :: itype

    ! The level where to add; 0=level independent part
    integer, intent(in) :: ilevel

!</input>

!<inputoutput>

    ! The collection
    type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<output>

    ! The pointer of the newly created value
    type(t_collctValue), pointer :: p_rvalue

!</output>
!</subroutine>

    ! local variables
    type(t_collctSection), pointer :: p_rSection
    type(t_collctLevel), pointer :: p_rLevel
    type(t_collctValue), dimension(:), pointer :: p_Rvalues
    integer :: i

    nullify(p_rvalue)

    ! The name must be given!
    if (sparamName .eq. '') return

    ! Get the section
    call collct_fetchsection(rcollection, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) return

    ! Get the level
    i = max(ilevel,0)
    if (i .gt. p_rsection%ilevelCount) return

    if (i .eq. 0) then
      p_rLevel => p_rsection%rlevel0
    else
      p_rLevel => p_rsection%p_Rlevels(i)
    end if

    ! Get space for the new structure in the list
    if (.not. associated(p_rLevel%p_Rvalues)) then

      ! Create an array for the entries
      allocate(p_rLevel%p_Rvalues(COLLCT_NVALUES))

      ! Take a pointer to the first position for storing data.
      p_rvalue => p_Rlevel%p_Rvalues(1)

    else if (p_rLevel%ivalueCount .ge. size(p_rLevel%p_Rvalues)) then

      ! Reallocate the entry array to get space for the new entry
      allocate(p_Rvalues(size(p_rLevel%p_Rvalues)+COLLCT_NVALUES))
      p_Rvalues(1:p_rLevel%ivalueCount) = p_rLevel%p_Rvalues (1:p_rLevel%ivalueCount)
      deallocate(p_rLevel%p_Rvalues)
      p_rLevel%p_Rvalues => p_Rvalues

      ! Store the value at the new free position.
      p_rvalue => p_Rlevel%p_Rvalues(p_rLevel%ivalueCount+1)

      ! There are definitely no holes anymore in the list on this
      ! level, as all free positions are filled up when we are here.
      p_rLevel%bisFull = .true.

    else
      ! Are there holes or can we directly take the maximum position?
      if (p_rLevel%bisFull) then
        ! No holes, take the next free position.
        i = p_rLevel%ivalueCount + 1
      else
        ! Find an empty position in the array.
        do i = 1, size(p_rLevel%p_Rvalues)
          if (p_rLevel%p_Rvalues(i)%itype .eq. COLLCT_UNDEFINED) exit
        end do

        ! This must work, otherwise ivalueCount is wrong!
        if (i .gt. size(p_rLevel%p_Rvalues)) then
          call output_line('Collection structure inconsistent!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'collct_addvalue')
          call sys_halt()
        end if

        ! If i is larger than ivalueCount, all holes are filled up again,
        ! i.e. we add 'behind' all elements.
        if (i .gt. p_rLevel%ivalueCount) p_rLevel%bisFull = .true.

      end if

      ! Take that value for storing data.
      p_rvalue => p_Rlevel%p_Rvalues(i)

    end if

    ! Fill the value with initial data and increase the counter.
    p_rLevel%ivalueCount = p_rLevel%ivalueCount + 1
    p_rvalue%itype = itype

    ! Store the name in upper case
    p_rvalue%sname = sys_upcase(sparamName)

  end subroutine collct_addvalue

  ! ***************************************************************************

!<subroutine>

  subroutine collct_initsection (rcollctSection, ssectionName)

!<description>

    ! Internal subroutine: Initialise a newly created section

!</description>

!<input>

    ! The name for the new section
    character(LEN=*), intent(in) :: ssectionName

!</input>

!<inputoutput>

    ! The section to initialise
    type(t_collctSection), intent(inout) :: rcollctSection

!</inputoutput>
!</subroutine>

    ! Simply allocate the pointers with an empty list
    allocate(rcollctSection%p_Rlevels(COLLCT_NVALUES))

    ! No variables in here for the moment
    rcollctSection%ilevelCount = 0

    ! Set the section name - uppercase
    rcollctSection%ssectionName = ssectionName
    call sys_toupper(rcollctSection%ssectionName)

    ! Init level 0
    call collct_initlevel (rcollctSection%rlevel0)

  end subroutine collct_initsection

  ! ***************************************************************************

!<subroutine>

  subroutine collct_donesection (rcollctSection)

!<description>

    ! Internal subroutine: Releases a section structure from memory.

!</description>

!<inputoutput>

    ! The section to release
    type(t_collctSection), intent(inout) :: rcollctSection

!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Clean up level 0
    call collct_donelevel (rcollctSection%rlevel0)

    ! Clean up the other levels in this section
    do i = rcollctSection%ilevelCount, 1, -1
      call collct_donelevel (rcollctSection%p_Rlevels(i))
    end do

    ! Deallocate the pointers with an empty list
    deallocate(rcollctSection%p_Rlevels)

    ! No variables in here for the moment
    rcollctSection%ilevelCount = 0

  end subroutine collct_donesection

  ! ***************************************************************************

!<subroutine>

  subroutine collct_reallocsection (rcollctSection, inewsize)

!<description>

    ! Internal subroutine: Reallocate a section.
    ! This increases the size of a section (i.e. the number of levels that
    ! can be stored in a section) by reallocation of the arrays.

!</description>

!<input>

    ! The new 'size' of the section, i.e. the new number of levels,
    ! the section should be able to handle.
    integer, intent(in) :: inewsize

!</input>

!<inputoutput>

    ! The section to reallocate.
    type(t_collctSection), intent(inout) :: rcollctSection

!</inputoutput>
!</subroutine>

    ! Pointers to new lists for replacing the old.
    type(t_collctLevel), dimension(:), pointer :: p_Rlevels

    ! local variables
    integer :: sz

    sz = max(size(rcollctSection%p_Rlevels),inewsize)

    if (size(rcollctSection%p_Rlevels) .eq. sz) return ! nothing to do

    ! Allocate the pointers for the new list
    allocate(p_Rlevels(sz))

    ! Copy the content of the old ones
    p_Rlevels(1:sz) = rcollctSection%p_Rlevels (1:sz)

    ! Throw away the old array, replace by the new one
    deallocate(rcollctSection%p_Rlevels)

    rcollctSection%p_Rlevels => p_Rlevels

  end subroutine collct_reallocsection

  ! ***************************************************************************

!<subroutine>

  subroutine collct_realloccollection (rcollection, inewsize)

!<description>

    ! Internal subroutine: Reallocate the collection.
    ! This increases the size of a a collection (i.e. the number of sections that
    ! can be stored in a collection) by reallocation of the arrays.

!</description>

!<input>

    ! The new 'size' of the section, i.e. the new number of parameters,
    ! the section should be able to handle.
    integer, intent(in) :: inewsize

!</input>

!<inputoutput>

    ! The section list to reallocate.
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine

    ! Pointers to new lists for replacing the old.
    type(t_collctSection), dimension(:), pointer :: p_Rsections

    ! local variables
    integer :: sz

    sz = max(size(rcollection%p_Rsections),inewsize)

    if (size(rcollection%p_Rsections) .eq. sz) return ! nothing to do

    ! Allocate the pointers for the new lists
    allocate(p_Rsections(sz))

    ! Copy the content of the old ones
    p_Rsections(1:sz) = rcollection%p_Rsections (1:sz)

    ! Throw away the old array, replace by the new one
    deallocate(rcollection%p_Rsections)

    rcollection%p_Rsections => p_Rsections

  end subroutine collct_realloccollection

  ! ***************************************************************************

!<subroutine>

  subroutine collct_fetchparameter_indir (rlevel, sname, iparamnum)

!<description>

    ! Internal subroutine: Search in a level subgroup for a parameter
    ! and return the index - or 0 if the parameter does not exist.

!</description>

!<input>

    ! The level where to search
    type(t_collctLevel), intent(in) :: rlevel

    ! The parameter name to look for.
    character(LEN=*), intent(in) :: sname

!</input>

!<output>

    ! The number of the parameter in the list or 0 if it does not exist.
    integer, intent(out) :: iparamnum

!</output>
!<subroutine>

    ! local variables
    integer :: i,nsections
    character(LEN=COLLCT_MLNAME)  :: sname2

    ! Convert the name to uppercase.
    sname2 = sname
    call sys_toupper (sname2)

    iparamnum = 0

    ! If the parameter list is empty, the section does not exist for sure
    if (rlevel%ivalueCount .eq. 0) return

    ! Loop through all sections to see if the section exists
    if (rlevel%bisfull) then
      nsections = rlevel%ivalueCount
    else
      nsections = size(rlevel%p_Rvalues)
    end if

    do i = 1, nsections
      if (rlevel%p_Rvalues(i)%sname .eq. sname2) then
        iparamnum = i
        return
      end if
    end do

  end subroutine collct_fetchparameter_indir

  ! ***************************************************************************

!<subroutine>

  subroutine collct_fetchparameter_direct (rcollection, ssectionName, ilevel, &
                                           sparameter, p_rvalue)

!<description>

    ! Internal subroutine: Search in a level subgroup for a parameter
    ! and return a pointer to it - or NULL if the parameter does not exist

!</description>

!<input>

    ! The collection where to search
    type(t_collection), intent(in) :: rcollection

    ! The section name where to search - '' identifies the unnamed section
    character(LEN=*), intent(in) :: ssectionName

    ! The level where to search; maybe 0 for level-independent parameters
    integer, intent(in) :: ilevel

    ! The parameter name to look for.
    character(LEN=*), intent(in) :: sparameter

!</input>

!<output>

    ! A pointer to the parameter.
    type(t_collctValue), pointer :: p_rvalue

!</output>
!</subroutine>

    ! local variables
    integer :: i,ilv
    type(t_collctSection), pointer :: p_rsection
    type(t_collctLevel), pointer :: p_rlevel

    ! Some basic checks
    if (rcollection%isectionCount .eq. 0) then
      call output_line('Collection not initalised!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_fetchparameter_direct')
      call sys_halt()
    end if

    nullify(p_rvalue)

    call collct_fetchsection(rcollection, ssectionName, p_rsection)

    if (.not. associated(p_rsection)) then
      ! Section does not exist - return NULL
      return
    end if

    ! Get the level - if it exists.
    ilv = max(0,ilevel)
    if(ilv .gt. p_rsection%ilevelCount) then
      ! Level does not exist - return NULL
      return
    end if

    if (ilv .eq. 0) then
      p_rlevel => p_rsection%rlevel0
    else
      p_rlevel => p_rsection%p_Rlevels(ilv)
    end if

    ! Get the index of the value
    call collct_fetchparameter_indir (p_rlevel, sparameter, i)

    ! and finally the pointer to it
    if (i .ne. 0) then
      p_rvalue => p_rlevel%p_Rvalues(i)
    end if

  end subroutine collct_fetchparameter_direct

  ! ***************************************************************************

!<subroutine>

  subroutine collct_fetchsection(rcollection, sname, p_rsection, isectionIndex)

!<description>

    ! Internal subroutine: Search in a collection for a section
    ! and return a pointer to the section - or NULL() if the section does
    ! not exist.

!</description>

!<input>

    ! The section.
    type(t_collection), intent(in) :: rcollection

    ! The parameter name to look for.
    character(LEN=*), intent(in) :: sname

!</input>

!<output>

    ! A pointer to the section or NULL() if it does not exist.
    type(t_collctSection), pointer :: p_rsection

    ! OPTIONAL: Index of the section in the section array.
    integer, intent(out), optional :: isectionIndex

!</output>
!</subroutine>

    ! local variables
    integer :: i
    character(LEN=COLLCT_MLSECTION) :: sname2

    ! Convert the name to uppercase.
    sname2 = sname
    call sys_toupper (sname2)

    nullify(p_rsection)

    ! If the parameter list is empty, the section does not exist for sure
    if (rcollection%isectionCount .eq. 0) return

    ! Loop through all sections to see if the section exists
    do i = 1, rcollection%isectionCount
      if (rcollection%p_Rsections(i)%ssectionName .eq. sname2) then
        p_rsection => rcollection%p_Rsections(i)

        if (present(isectionIndex)) isectionIndex = i

        return
      end if
    end do

  end subroutine collct_fetchsection

  ! ***************************************************************************

!<subroutine>

  subroutine collct_fetchlevel(rcollection, ssectionName, ilevel, p_rlevel)

!<description>

    ! Internal subroutine: Search in a collection for a level in a section
    ! and return a pointer to the level - or NULL)( if the level or the section
    ! does not exist).

!</description>

!<input>

    ! The section.
    type(t_collection), intent(in) :: rcollection

    ! The section name to look for.
    character(LEN=*), intent(in) :: ssectionName

    ! The level number where to search for
    integer, intent(in) :: ilevel

!</input>

!<output>

    ! A pointer to the level - or NULL() if the level does not exist
    type(t_collctLevel), pointer :: p_rlevel

!</output>
!</subroutine>

    ! local variables
    integer :: i
    type(t_collctSection), pointer :: p_rsection

    nullify(p_rlevel)

    ! Get the section
    call collct_fetchsection(rcollection, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) return

    ! Get the level
    i = max(0,ilevel)
    if (i .gt. p_rsection%ilevelCount) return

    if (i .eq. 0) then
      ! Get level 0 - the always-existing level
      p_rlevel => p_rsection%rlevel0
    else
      p_rlevel => p_rsection%p_Rlevels(i)
    end if

  end subroutine collct_fetchlevel

  ! ***************************************************************************

!<subroutine>

  subroutine collct_init (rcollection)

!<description>

  ! This routine initialises a collection. It must be applied
  ! before doing anything to it, just to initialise.

!</description>

!<inputoutput>

  ! The parameter list to initialise.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!</subroutine>

    ! Allocate memory for the sections. We have at least one unnamed section
    ! and allocate memory for COLLCT_NSECTIONS more named sections in advance.
    allocate(rcollection%p_Rsections(1+COLLCT_NSECTIONS))

    ! Add one section - the unnamed one
    rcollection%isectionCount = 1
    call collct_initsection (rcollection%p_Rsections(1), '')

  end subroutine collct_init

  ! ***************************************************************************

!<subroutine>

  subroutine collct_done (rcollection)

!<description>

  ! This routine releases a parameter list. All memory allocated by the
  ! parameter list is released.

!</description>

!<inputoutput>

  ! The parameter list to release.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Structure initialised?
    if (rcollection%isectionCount .eq. 0) then
      call output_line('Trying to release an uninitialised collection.',&
                       OU_CLASS_WARNING,OU_MODE_STD,'collct_done')
      return
    end if

    ! Loop through the sections and release the content
    do i = rcollection%isectionCount, 1, -1
      call collct_donesection(rcollection%p_Rsections(i))
    end do

    ! Release all section structures
    deallocate(rcollection%p_Rsections)

    ! Mark the structure as 'empty', finish
    rcollection%isectionCount = 0

  end subroutine collct_done

  ! ***************************************************************************

!<function>

  integer function collct_getmaxlevel_direct (rcollection,ssectionName) result (ilevel)

!<description>
  ! Returns the maximum available level in section ssectionname of a collection.
!</description>

!<input>

  ! The collection
  type(t_collection), intent(in) :: rcollection

  ! The name of the section
  character(LEN=*), intent(in) :: ssectionName

!</input>

!<result>
  ! The maximum level in rcollection.
!</result>

!</function>

    type(t_collctSection), pointer :: p_rsection

    ! Fetch the section number
    call collct_fetchsection(rcollection, ssectionName, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line('Section not found: '//ssectionName,&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_getmaxlevel_direct')
      call sys_halt()
    end if

    ilevel = p_rsection%ilevelCount

  end function collct_getmaxlevel_direct

  ! ***************************************************************************

!<function>

  integer function collct_getmaxlevel_indir (rsection) result (ilevel)

!<description>
  ! Returns the maximum available level in section rsection.
  ! A return value of 0 indicates that only level-independent information
  ! is available in the block at the moment.
!</description>

!<input>

  ! The collection
  type(t_collctsection), intent(in) :: rsection

!</input>

!<result>
  ! The maximum level in rcollection.
!</result>

!</function>

    ilevel = rsection%ilevelCount

  end function collct_getmaxlevel_indir

  ! ***************************************************************************

!<subroutine>

  subroutine collct_addlevel_indir (rsection, ilevelid)

!<description>
  ! Adds a new level to a section rsection.
  ! ilevelid returns the id of the new level and is >= 1.
!</description>

!<inputoutput>

  ! The section where to add the level
  type(t_collctSection), intent(inout) :: rsection

!</inputoutput>

!<output>

  ! OPTIONAL: The id of the newly added level.
  integer, intent(out), optional :: ilevelid

!</output>

!</subroutine>

    ! Add a new level - reallocate the level list if necessary
    if (.not. associated(rsection%p_Rlevels)) then
      allocate(rsection%p_Rlevels(COLLCT_NLEVELS))
    else if (rsection%ilevelCount .eq. size(rsection%p_Rlevels)) then
      call collct_reallocsection (rsection, size(rsection%p_Rlevels)+COLLCT_NLEVELS)
    end if
    rsection%ilevelCount = rsection%ilevelCount + 1

    ! Initialise the new level.
    call collct_initlevel(rsection%p_Rlevels(rsection%ilevelCount))

    if (present(ilevelid)) ilevelid = rsection%ilevelCount

  end subroutine collct_addlevel_indir

  ! ***************************************************************************

!<subroutine>

  subroutine collct_addlevel_direct (rcollection, ilevelid, ssectionName)

!<description>
  ! Adds a new level to a section with the name ssectionName.
  ! If ssectionName='', the new level is added to the unnamed section.
  ! ilevelid returns the id of the new level and is >= 1.
!</description>

!<inputoutput>

  ! The collection where to add the level
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! OPTIONAL: The section name where to add the level. If ='' or not
  ! given, the level is added to the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! The id of the newly added level.
  integer, intent(out) :: ilevelid

!</output>

!</subroutine>

    ! local variables
    type(t_collctSection), pointer :: p_rsection

    if (rcollection%isectionCount .eq. 0) then
      call output_line('Collection not initalised!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_addlevel_direct')
      call sys_halt()
    end if

    ! Get the section
    if (present(ssectionName)) then
      call collct_fetchsection(rcollection, ssectionName, p_rsection)
      if (.not. associated(p_rsection)) then
        call output_line('Section not found: '//ssectionName,&
                         OU_CLASS_ERROR,OU_MODE_STD,'collct_addlevel_direct')
        call sys_halt()
      end if
    else
      ! unnamed section
      p_rsection => rcollection%p_Rsections(1)
    end if

    ! Add the level
    call collct_addlevel_indir (p_rsection, ilevelid)

  end subroutine collct_addlevel_direct

  ! ***************************************************************************

!<subroutine>

  subroutine collct_addlevel_all (rcollection)

!<description>
  ! Adds a new level to all sections in the collection.
!</description>

!<inputoutput>

  ! The collection where to add the level
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,ilevelid

    if (rcollection%isectionCount .eq. 0) then
      call output_line('Collection not initalised!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_addlevel_all')
      call sys_halt()
    end if

    ! Go through all sections
    do i = 1, rcollection%isectionCount
      ! Add the level there
      call collct_addlevel_indir (rcollection%p_Rsections(i), ilevelid)
    end do

  end subroutine collct_addlevel_all

  ! ***************************************************************************

!<subroutine>

  subroutine collct_addsection (rcollection, ssectionName)

!<description>
  ! Adds a new section to a collection. The section gets the name ssectionName.
!</description>

!<inputoutput>

  ! The parameter list where to add the section.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The name of the section. Must not be ''!
  character(LEN=*), intent(in) :: ssectionName

!</input>

!</subroutine>

    if (ssectionName .eq. '') then
      call output_line('Empty section name!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_addsection')
      call sys_halt()
    end if

    if (rcollection%isectionCount .eq. 0) then
      call output_line('Collection not initalised!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_addsection')
      call sys_halt()
    end if

    ! Add a new section - reallocate the memory if necessary
    if (.not. associated(rcollection%p_Rsections)) then
      allocate(rcollection%p_Rsections(COLLCT_NSECTIONS))
    else if (rcollection%isectionCount .eq. size(rcollection%p_Rsections)) then
      call collct_realloccollection (rcollection, size(rcollection%p_Rsections)+COLLCT_NSECTIONS)
    end if
    rcollection%isectionCount = rcollection%isectionCount + 1

    ! Initialise the new section.
    call collct_initsection(rcollection%p_Rsections(rcollection%isectionCount),ssectionName)

  end subroutine collct_addsection

  ! ***************************************************************************

!<subroutine>

  subroutine collct_deletesection (rcollection, ssectionName)

!<description>
  ! Removes the section ssectionName with all its content from the collection.
  !
  ! Warning: Removing a section from the collection makes all pointers to
  ! sections invalid!
!</description>

!<inputoutput>

  ! The parameter list where to add the section.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The name of the section. Must not be ''!
  character(LEN=*), intent(in) :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctSection), pointer :: p_rsection
    integer :: isectionIndex,i

    ! Some basic checks
    if (rcollection%isectionCount .eq. 0) then
      call output_line('Collection not initalised!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_deletesection')
      call sys_halt()
    end if

    call collct_fetchsection(rcollection, ssectionName, p_rsection, isectionIndex)

    if (.not. associated(p_rsection)) then
      call output_line('Section does not exist: '//ssectionName,&
                       OU_CLASS_WARNING,OU_MODE_STD,'collct_deletesection')
      return
    end if

    ! Remove the section data
    call collct_donesection (p_rsection)

    ! Now 'relocate' all sections. This makes all pointers to sections
    ! in the main application invalid, but it is rather unlikely that the
    ! application maintains a section pointer in a stage where sections are
    ! modified - hopefully :-)

    do i = isectionIndex+1, rcollection%isectionCount
      rcollection%p_Rsections(i-1) = rcollection%p_Rsections(i)
    end do

    ! Reduce the number of sections in the collection, finish.
    rcollection%isectionCount = rcollection%isectionCount - 1

  end subroutine collct_deletesection

  ! ***************************************************************************

!<subroutine>

  subroutine collct_cleanupvalue (rvalue)

!<description>
  ! Cleans up a value structure. Releases all associated memory.
!</description>

!<inputoutput>
  ! The value to clean up.
  type(t_collctValue), intent(inout) :: rvalue
!</inputoutput>

!</subroutine>

    if (rvalue%itype .ne. COLLCT_UNDEFINED) then
      ! Deleting it is simple: Clear the name and assign all values
      ! to default values. Then the node is 'dead' and can be reused
      ! later.

      rvalue%sname = ''
      rvalue%csvalue = " "
      if (associated(rvalue%p_svalue)) deallocate(rvalue%p_svalue)
      rvalue%ivalue = 0
      rvalue%dvalue = 0.0_DP
      nullify(rvalue%p_rdiscretisation)
      nullify(rvalue%p_rtriangulation)
      nullify(rvalue%p_rvectorScalar)
      nullify(rvalue%p_rmatrixScalar)
      nullify(rvalue%p_rvector)
      nullify(rvalue%p_rmatrix)
      nullify(rvalue%p_rparlist)
      nullify(rvalue%p_rlinearSolver)
      nullify(rvalue%p_rdiscreteBC)
      nullify(rvalue%p_rboundary)
      nullify(rvalue%p_rboundaryConditions)
      nullify(rvalue%p_rparser)
      nullify(rvalue%p_RfilterChain)
      nullify(rvalue%p_rcollection)
      nullify(rvalue%p_rhadapt)
      nullify(rvalue%p_rafcstab)

      if (associated(rvalue%p_Iarray)) deallocate(rvalue%p_Iarray)
      if (associated(rvalue%p_Darray)) deallocate(rvalue%p_Darray)

      rvalue%itype = COLLCT_UNDEFINED
    end if

  end subroutine collct_cleanupvalue

  ! ***************************************************************************

!<subroutine>

  subroutine collct_deletevalue (rcollection, sparameter, ilevel, ssectionName)

!<description>
  ! Deletes the stored parameter 'sparameter' from the collection.
  ! This simply removes a pointer or value from the collection, no memory
  ! is released.
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>
  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName
!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    type(t_collctLevel), pointer :: p_rlevel
    integer :: ilv

    ilv = 0
    if (present(ilevel)) ilv = ilevel

    ! Get the parameter
    if (present(ssectionName)) then
      call collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                         sparameter, p_rvalue)
    else
      call collct_fetchparameter_direct (rcollection, '', ilv, &
                                         sparameter, p_rvalue)
    end if

    ! Does the value exist?
    if (associated(p_rvalue)) then

      ! Ok, we have the parameter. Clean it up.
      call collct_cleanupvalue(p_rvalue);

      ! Modify the level info:
      if (present(ssectionName)) then
        call collct_fetchlevel(rcollection, ssectionName, ilv, p_rlevel)
      else
        call collct_fetchlevel(rcollection, '', ilv, p_rlevel)
      end if

      ! Decrement the value counter
      p_rLevel%ivalueCount = p_rLevel%ivalueCount - 1

      ! Indicate that there is at least one free position now
      p_rLevel%bisFull = .false.

    end if

  end subroutine collct_deletevalue

  ! ***************************************************************************

!<subroutine>

  subroutine collct_printStatistics (rcollection)

!<description>
  ! This routine prints the current content of the whole collection in
  ! user-readable form to screen.
!</description>

!<input>
  ! The parameter list.
  type(t_collection), intent(in) :: rcollection
!</input>

!</subroutine>

    ! local variables
    integer :: isection, ilevel

    ! Loop through all sections in the collection
    do isection = 1, rcollection%isectionCount

      ! Print the segment name:
      call output_line ('['//trim(rcollection%p_Rsections(isection)%ssectionName)//']')

      ! Print the content of level 0 of the current section
      call printlevel(0,rcollection%p_Rsections(isection)%rlevel0)

      ! Loop through all levels and print them
      do ilevel = 1,rcollection%p_Rsections(isection)%ilevelCount
        call printlevel(ilevel,rcollection%p_Rsections(isection)%p_Rlevels(ilevel))
      end do

    end do

  contains

    !****************************************************************
    ! Print out the content of a level to screen
    subroutine printlevel (ilevel,rlevel)

      ! Number of the level; 0=level 0
      integer, intent(in) :: ilevel

      ! The level structure
      type(t_collctLevel), intent(in) :: rlevel

      ! local variables
      integer :: ivalue
      type(t_collctValue), pointer :: p_rvalue

      ! Is there data in the level?
      if (rlevel%ivalueCount .ne. 0) then

        ! A level > 0 has a name:
        if (ilevel .gt. 0) then
          call output_line ('  [Level '//trim(sys_siL(ilevel,8))//']')
        end if

        ! Loop through all values
        do ivalue = 1,size(rlevel%p_Rvalues)
          ! Get the value structure
          p_rvalue => rlevel%p_Rvalues(ivalue)

          ! Is there data in it?
          if (p_rvalue%itype .ne. COLLCT_UNDEFINED) then

            ! Print that there is data:
            call output_line ('  '//trim(p_rvalue%sname)//' (Type: '//&
                              trim(sys_siL(p_rvalue%itype,8))//')')

          end if
        end do

      end if

    end subroutine printlevel

  end subroutine collct_printStatistics

  ! ***************************************************************************

!<function>

  integer function collct_queryvalue_indir (rsection, sparameter, ilevel) &
               result (exists)

!<description>
  ! Checks whether a parameter sparameter exists in in section rsection.
!</description>

!<result>
  ! >0, if the parameter exists in the section rsection (it is the index of the
  !     parameter within the section),
  ! =0, if the parameter does not exist whithin the section.
!</result>

!<input>

  ! The parameter list.
  type(t_collctSection), intent(in), target :: rsection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

!</input>

!</function>

    ! local variables
    integer :: ilv
    type(t_collctLevel), pointer :: p_rlevel

    exists = 0

    if (.not. present(ilevel)) then
      ilv = 0
    else
      ilv = max(ilevel,0)
    end if

    ! Get the right level structure
    if (ilv .eq. 0) then
      p_rlevel => rsection%rlevel0
    else if (ilevel .gt. rsection%ilevelCount) then
      call output_line('Level out of bounds!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_queryvalue_indir')
      call sys_halt()
    else
      p_rlevel => rsection%p_Rlevels(ilv)
    end if

    ! Search for the parameter in the level
    call collct_fetchparameter_indir(p_rlevel, sparameter, exists)

  end function collct_queryvalue_indir

  ! ***************************************************************************

!<function>

  integer function collct_queryvalue_direct (rcollection, sparameter, &
                        ilevel, ssectionName) result (exists)

!<description>
  ! Checks whether a parameter sparameter exists in in section.
  ! ilevel and ssectionName are optional parameters.
  ! If ilevel=0 or not given, the parameter is searched in the level-independent
  ! parameter set, otherwise on level ilevel.
  ! If ssectionName='' or not given, the parameter is searched in the unnamed
  ! section, otherwise in section ssectionName.
!</description>

!<result>
  ! >0, if the parameter exists in the section (it is the index of the
  !     parameter within the section),
  ! =0, if the parameter does not exist whithin the section.
!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(in) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</function>

    ! local variables
    integer :: ilv
    type(t_collctSection), pointer :: p_rsection

    if (rcollection%isectionCount .eq. 0) then
      call output_line('Collection not initalised!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'collct_queryvalue_direct')
      call sys_halt()
    end if

    ! Fall back to the 'indirect' search.
    if (.not. present(ilevel)) then
      ilv = 0
    else
      ilv = max(0,ilevel)
    end if

    if (.not. present(ssectionName)) then
      p_rsection => rcollection%p_Rsections(1)
    else
      call collct_fetchsection(rcollection, ssectionName, p_rsection)
      if (.not. associated(p_rsection)) then
        call output_line('Section not found',&
                         OU_CLASS_ERROR,OU_MODE_STD,'collct_queryvalue_direct')
        call sys_halt()
      end if
    end if

    exists = collct_queryvalue_indir (p_rsection, sparameter, ilv)

  end function collct_queryvalue_direct

  ! ***************************************************************************

!<function>

  integer function collct_gettype (rcollection, sparameter, &
                        ilevel, ssectionName, bexists) result (itype)
!<description>
  ! Returns the data type of the parameter sparameter.
!</description>

!<result>
  ! One of the COLLCT_XXXX constants (COLLCT_UNDEFINED,...).
!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(in) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    integer :: ilv

    ilv = 0
    if (present(ilevel)) ilv = ilevel

    ! Get the parameter
    if (present(ssectionName)) then
      call collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                         sparameter, p_rvalue)
    else
      call collct_fetchparameter_direct (rcollection, '', ilv, &
                                         sparameter, p_rvalue)
    end if

    ! Return whether or not that thing exists
    if (present(bexists)) bexists = associated(p_rvalue)

    ! Return the quantity
    if (associated(p_rvalue)) then
      itype = p_rvalue%itype
    else
      itype = COLLCT_UNDEFINED
    end if

  end function collct_gettype

  ! ***************************************************************************

!<function>

  integer function collct_gettag (rcollection, sparameter, &
                        ilevel, ssectionName, bexists) result (itag)
!<description>
  ! Returns the user defined tag of the parameter sparameter.
!</description>

!<result>
  ! User defined integer tag.
!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(in) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    integer :: ilv

    ilv = 0
    if (present(ilevel)) ilv = ilevel

    ! Get the parameter
    if (present(ssectionName)) then
      call collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                         sparameter, p_rvalue)
    else
      call collct_fetchparameter_direct (rcollection, '', ilv, &
                                         sparameter, p_rvalue)
    end if

    ! Return whether or not that thing exists
    if (present(bexists)) bexists = associated(p_rvalue)

    ! Return the quantity
    if (associated(p_rvalue)) then
      itag = p_rvalue%itag
    else
      itag = 0
    end if

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine collct_settag (rcollection, sparameter, itag, &
                        ilevel, ssectionName, bexists)
!<description>
  ! Sets the user defined tag of the parameter sparameter.
!</description>

!<result>
  ! User defined integer tag.
!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(in) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! Tag that should be associated to the parameter.
  ! User defined integer value.
  integer, intent(in) :: itag

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    integer :: ilv

    ilv = 0
    if (present(ilevel)) ilv = ilevel

    ! Get the parameter
    if (present(ssectionName)) then
      call collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                         sparameter, p_rvalue)
    else
      call collct_fetchparameter_direct (rcollection, '', ilv, &
                                         sparameter, p_rvalue)
    end if

    ! Return whether or not that thing exists
    if (present(bexists)) bexists = associated(p_rvalue)

    ! Return the quantity
    if (associated(p_rvalue)) then
      p_rvalue%itag = itag
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine collct_getvalue_struc (rcollection, sparameter, itype, &
                                    badd, p_rvalue, ilevel, bexists, ssectionName)
!<description>
  ! INTERNAL SUBROUTINE!
  ! Searches in the collection rcollection on the level ilevel for the
  ! parameter sparameter. Returns a pointer to the value structure in
  ! p_rvalue.
  !
  ! If the value does not exists and badd=false, an error is thrown.
  ! If the value does not exists and badd=true, a new value is created.
  ! However, no memory is allocated for array-type values.
  !
  ! If bexists is present, it is set to TRUE/FALSE to indicate whether
  ! the value exists at all. If bexists does not exist, an error is thrown.
  ! ssectionName allows to specify the name of a section where to search
  ! for the value structure.
  !
  ! itype specifies the desired type of the value.
  ! If itype <> COLLCT_UNDEFINED is specified, there is a type check:
  ! If the value exists but is not of type itype, an error is thrown.
!</description>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! Expected type for that parameter.
  ! =COLLCT_UNDEFINED: suppress type check.
  integer, intent(in) :: itype

  ! Whether to add the parameter if it does not exists.
  ! =TRUE: add the parameter if necessary.
  ! =FALSE: throw an error that the parameter does not exist.
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

  ! A pointer to the structure with name sparameter.
  ! Points to NULL() if the value does not exist.
  type(t_collctValue), pointer :: p_rvalue

!</output>

!</subroutine>

    ! local variables
    integer :: ilv

    ilv = 0
    if (present(ilevel)) ilv = ilevel

    ! Get the parameter
    if (present(ssectionName)) then
      call collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                         sparameter, p_rvalue)
    else
      call collct_fetchparameter_direct (rcollection, '', ilv, &
                                         sparameter, p_rvalue)
    end if

    ! Return whether or not that thing exists.
    ! If it does not exist, create it or throw an error.
    if (present(bexists)) then
      bexists = associated(p_rvalue)
    else
      if ((.not. badd) .and. (.not. associated(p_rvalue))) then
        call output_line('Parameter '//trim(sparameter)//' not found!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'collct_getvalue_struc')
        call sys_halt()
      end if
    end if

    if (associated(p_rvalue)) then

      ! Type check
      if (itype .ne. COLLCT_UNDEFINED) then
        ! Throw an error if the type is wrong. Otherwise, get the value.
        if (p_rvalue%itype .ne. itype) then
          call output_line('Parameter: '//trim(sparameter)//&
                           ' is of type '//trim(sys_siL(p_rvalue%itype,5))//&
                           ' but expected as '//trim(sys_siL(itype,5))//'!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'collct_getvalue_struc')
          call sys_halt()
        end if
      end if

    else
      if (badd) then
        if (present(ssectionName)) then
          call collct_addvalue (rcollection, ssectionName, sparameter, &
                                itype, ilv, p_rvalue)
        else
          call collct_addvalue (rcollection, '', sparameter, &
                                itype, ilv, p_rvalue)
        end if
        if (present(bexists)) bexists = .true.
      end if
    end if

  end subroutine collct_getvalue_struc

  ! ***************************************************************************

!<function>

  character function collct_getvalue_char (rcollection, sparameter, &
                                           ilevel, ssectionName, bexists) result (value)
!<description>
  ! Returns the the parameter sparameter as character.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_CHARACTER,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value = p_rvalue%csvalue
    else
      value = " "
    end if

  end function collct_getvalue_char

  ! ***************************************************************************

!<subroutine>

  subroutine collct_getvalue_string (rcollection, sparameter, value, &
                                     ilevel, ssectionName, bexists)
!<description>
  ! Returns the the parameter sparameter as string.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! The string.
  character(LEN=*) :: value

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    integer :: i,j

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_STRING,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      i = min(size(p_rvalue%p_svalue),len(value))
      value = ''
      do j=1,i
        value(j:j) = p_rvalue%p_svalue(j)
      end do
    else
      value = ''
    end if

  end subroutine collct_getvalue_string

  ! ***************************************************************************

!<function>

  integer function collct_getvalue_int (rcollection, sparameter, &
                                        ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as integer.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGER,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value = p_rvalue%ivalue
    else
      value = 0
    end if

  end function collct_getvalue_int

  ! ***************************************************************************

!<subroutine>

  subroutine collct_getvalue_intarr (rcollection, sparameter, value, &
                                     ilevel, ssectionName, bexists)
!<description>
  ! Returns the the parameter sparameter as integer array.
  ! An error is thrown if the value is of the wrong type.
  !
  ! The routine returns the actual array, not a reference to it. The previous
  ! content of the destination array is overwritten.
!</description>

!<inputoutput>
  ! The value of the parameter an array. The destination array must have the
  ! same length as the array in the collection!
  ! A standard value if the value does not exist.
  integer, dimension(:), intent(inout) :: value
!</inputoutput>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>
  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists
!</output>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGERARR,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      if (size(p_rvalue%p_Iarray) .ne. size(value)) then
        call output_line ('Destination array has the wrong length!',&
                          OU_CLASS_ERROR,OU_MODE_STD,'collct_getvalue_intarr')
        call sys_halt()
      end if
      value = p_rvalue%p_Iarray
    else
      value = 0
    end if

  end subroutine collct_getvalue_intarr

  ! ***************************************************************************

!<subroutine>

  subroutine collct_getvalue_iarrp (rcollection, sparameter, value, &
                                    ilevel, ssectionName, bexists)
!<description>
  ! Returns the the parameter sparameter as integer array.
  ! An error is thrown if the value is of the wrong type.
  !
  ! The routine returns a reference to the stored array.
!</description>

!<inputoutput>
  ! The value of the parameter an array. The destination array must have the
  ! same length as the array in the collection!
  ! A standard value if the value does not exist.
  integer, dimension(:), pointer :: value
!</inputoutput>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>
  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists
!</output>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGERARRP,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    value => p_rvalue%p_Iarray

  end subroutine

  ! ***************************************************************************

!<function>

  real(DP) function collct_getvalue_real (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as real.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_REAL,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value = p_rvalue%dvalue
    else
      value = 0
    end if

  end function collct_getvalue_real

  ! ***************************************************************************

!<subroutine>

  subroutine collct_getvalue_realarr (rcollection, sparameter, value, &
                                  ilevel, ssectionName, bexists)
!<description>
  ! Returns the the parameter sparameter as real array.
  ! An error is thrown if the value is of the wrong type.
  !
  ! The routine returns the actual array, not a reference to it. The previous
  ! content of the destination array is overwritten.
!</description>

!<inputoutput>
  ! The value of the parameter an array. The destination array must have the
  ! same length as the array in the collection!
  ! A standard value if the value does not exist.
  real(DP), dimension(:), intent(inout) :: value
!</inputoutput>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_REALARR,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      if (size(p_rvalue%p_Darray) .ne. size(value)) then
        call output_line ('Destination array has the wrong length!',&
                          OU_CLASS_ERROR,OU_MODE_STD,'collct_getvalue_realarr')
        call sys_halt()
      end if
      value = p_rvalue%p_Darray
    else
      value = 0
    end if

  end subroutine collct_getvalue_realarr

  ! ***************************************************************************

!<subroutine>

  subroutine collct_getvalue_rarrp (rcollection, sparameter, value, &
                                  ilevel, ssectionName, bexists)
!<description>
  ! Returns the the parameter sparameter as real array.
  ! An error is thrown if the value is of the wrong type.
  !
  ! The routine returns a reference to the array.
!</description>

!<inputoutput>
  ! The value of the parameter an array. The destination array must have the
  ! same length as the array in the collection!
  ! A standard value if the value does not exist.
  real(DP), dimension(:), pointer :: value
!</inputoutput>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_REALARRP,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    value => p_rvalue%p_Darray

  end subroutine

  ! ***************************************************************************

!<function>

  function collct_getvalue_discr (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a discretisation
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>

  type(t_spatialDiscretisation), pointer :: value

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCR,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rdiscretisation
    else
      nullify(value)
    end if

  end function collct_getvalue_discr

  ! ***************************************************************************

!<function>

  function collct_getvalue_bldiscr (rcollection, sparameter, &
                                    ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a block discretisation
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_blockDiscretisation), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BLDISCR,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rblockDiscretisation
    else
      nullify(value)
    end if

  end function collct_getvalue_bldiscr

  ! ***************************************************************************

!<function>

  function collct_getvalue_tria (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a triangulation
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_triangulation), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_TRIA,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rtriangulation
    else
      nullify(value)
    end if

  end function collct_getvalue_tria

  ! ***************************************************************************

!<function>

  function collct_getvalue_bdry (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a boundary
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
    type(t_boundary), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARY,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rboundary
    else
      nullify(value)
    end if

  end function collct_getvalue_bdry

  ! ***************************************************************************

!<function>

  function collct_getvalue_bdreg (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a boundary region
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
    type(t_boundaryRegion), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BDREGION,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rboundaryRegion
    else
      nullify(value)
    end if

  end function collct_getvalue_bdreg

  ! ***************************************************************************

!<function>

  function collct_getvalue_bc (rcollection, sparameter, &
                               ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a triangulation
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_boundaryConditions), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARYCOND,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rboundaryConditions
    else
      nullify(value)
    end if

  end function collct_getvalue_bc

  ! ***************************************************************************

!<function>

  function collct_getvalue_vecsca (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar vector.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_vectorScalar), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAVECTOR,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rvectorScalar
    else
      nullify(value)
    end if

  end function collct_getvalue_vecsca

  ! ***************************************************************************

!<function>

  function collct_getvalue_matsca (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar matrix.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_matrixScalar), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAMATRIX,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rmatrixScalar
    else
      nullify(value)
    end if

  end function collct_getvalue_matsca

  ! ***************************************************************************

!<function>

  function collct_getvalue_vec (rcollection, sparameter, &
                                ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a block vector.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_vectorBlock), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKVECTOR,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rvector
    else
      nullify(value)
    end if

  end function collct_getvalue_vec

  ! ***************************************************************************

!<function>

  function collct_getvalue_mat (rcollection, sparameter, &
                               ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar matrix.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_matrixBlock), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKMATRIX,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rmatrix
    else
      nullify(value)
    end if

  end function collct_getvalue_mat

  ! ***************************************************************************

!<function>

  function collct_getvalue_linsol (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a linear solver
  ! configuration node.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_linsolNode), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_LINSOL,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rlinearSolver
    else
      nullify(value)
    end if

  end function collct_getvalue_linsol

  ! ***************************************************************************

!<function>

  function collct_getvalue_discbc (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a discrete boundary
  ! condition structure
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_discreteBC), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCRBC,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rdiscreteBC
    else
      nullify(value)
    end if

  end function collct_getvalue_discbc

  ! ***************************************************************************

!<function>

  function collct_getvalue_parlst (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a parameter list.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_parlist), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_PARAMETERS,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rparlist
    else
      nullify(value)
    end if

  end function collct_getvalue_parlst

  ! ***************************************************************************

!<function>

  function collct_getvalue_ilvpsc (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar interlevel
  ! projection structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_interlevelProjectionScalar), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJSC,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rilvProjectionSc
    else
      nullify(value)
    end if

  end function collct_getvalue_ilvpsc

  ! ***************************************************************************

!<function>

  function collct_getvalue_ilvp (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a (block) interlevel
  ! projection structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_interlevelProjectionBlock), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJ,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rilvProjection
    else
      nullify(value)
    end if

  end function collct_getvalue_ilvp

  ! ***************************************************************************

!<function>

  function collct_getvalue_coll (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a collection object.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_collection), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_COLLECTION,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rcollection
    else
      nullify(value)
    end if

  end function collct_getvalue_coll

  ! ***************************************************************************

!<function>

  function collct_getvalue_pars (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a parser object.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_fparser), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FPARSER,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rparser
    else
      nullify(value)
    end if

  end function collct_getvalue_pars

  ! ***************************************************************************

!<function>

  function collct_getvalue_geom (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a geometry object.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_geometryObject), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GEOMETRY,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rgeometry
    else
      nullify(value)
    end if

  end function collct_getvalue_geom

  ! ***************************************************************************

!<function>

  function collct_getvalue_particles(rcollection, sparameter, &
                                     ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a particle collection object.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_particleCollection), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_PARTICLES,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rparticles
    else
      nullify(value)
    end if

  end function collct_getvalue_particles

  ! ***************************************************************************

!<function>

  function collct_getvalue_particles3D(rcollection, sparameter, &
                                       ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a particle3D collection object.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_particleCollection3D), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_PARTICLES3D,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rparticles3D
    else
      nullify(value)
    end if

  end function collct_getvalue_particles3D

  ! ***************************************************************************

!<function>

  function collct_getvalue_fchn (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a filter chain.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_filterChain), dimension(:), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FILTERCHAIN,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_RfilterChain
    else
      nullify(value)
    end if

  end function collct_getvalue_fchn

  ! ***************************************************************************

!<function>

  function collct_getvalue_hadapt (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a grid adaptation structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_hadapt), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_HADAPT,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rhadapt
    else
      nullify(value)
    end if

  end function collct_getvalue_hadapt

  ! ***************************************************************************

!<function>

  function collct_getvalue_afcstab (rcollection, sparameter, &
                                    ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a stabilisation structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_afcstab), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_AFCSTAB,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rafcstab
    else
      nullify(value)
    end if

  end function collct_getvalue_afcstab

  ! ***************************************************************************

!<function>

  function collct_getvalue_timer (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a timer structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_timer), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_TIMER,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rtimer
    else
      nullify(value)
    end if

  end function collct_getvalue_timer

  ! ***************************************************************************

!<function>

  function collct_getvalue_graph (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a graph structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_graph), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GRAPH,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rgraph
    else
      nullify(value)
    end if

  end function collct_getvalue_graph

  ! ***************************************************************************

!<function>

  function collct_getvalue_mshh (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a mesh hierarchy.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_meshHierarchy), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_MSHHIERARCHY,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rmeshhierarchy
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_fesp (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a finite element space level.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_feSpaceLevel), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FESPACE,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rfeSpace
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_feh (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a finite element hierarchy.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_feHierarchy), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FEHIERARCHY,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rfeHierarchy
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_tsh (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a time scale hierarchy.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_timescaleHierarchy), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_TSHIERARCHY,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rtsHierarchy
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_mlprjh (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a multilevel projection structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_interlevelProjectionHier), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_MLPRJHIERARCHY,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rmlprjHierarchy
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_gfemset (rcollection, sparameter, &
                                    ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a group finite element set.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_groupFEMSet), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GROUPFEMSET,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rgroupFEMSet
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_gfemblk (rcollection, sparameter, &
                                    ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a group finite element block.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_groupFEMBlock), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GROUPFEMBLOCK,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rgroupFEMBlock
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_cubinfo (rcollection, sparameter, &
                                    ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar cubature info structure.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_scalarCubatureInfo), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_CUBINFO,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rcubatureInfo
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<function>

  function collct_getvalue_genobj (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) result(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a generic object.
  ! An error is thrown if the value is of the wrong type.
!</description>

!<result>

  ! The value of the parameter.
  ! A standard value if the value does not exist.
  type(t_genericObject), pointer :: value

!</result>

!<input>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There is no error thrown if a variable does not exist.
  logical, intent(out), optional :: bexists

!</output>

!</function>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue

    ! Get the pointer to the parameter
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GENOBJ,&
                                .false., p_rvalue, ilevel, bexists, ssectionName)

    ! Return the quantity
    if (associated(p_rvalue)) then
      value => p_rvalue%p_rgenericObject
    else
      nullify(value)
    end if

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_char (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  character, intent(in) :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_CHARACTER,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%csvalue = value

  end subroutine collct_setvalue_char

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_string (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  character(LEN=*), intent(in) :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists
    integer :: i

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_STRING,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value; slightly more complicated for strings
    allocate(p_rvalue%p_svalue(len(value)))
    do i = 1, len(value)
      p_rvalue%p_svalue(i) = value (i:i)
    end do

  end subroutine collct_setvalue_string

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_int (rcollection, sparameter, value, badd, &
                                  ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  integer, intent(in) :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGER,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%ivalue = value

  end subroutine collct_setvalue_int

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_intarr (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
  !
  ! The routine creates a new memory block on the heap and copies the
  ! data to there. There is no reference to the old array stored.
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  integer, dimension(:), intent(in) :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGERARR,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! (Re-)allocate memory for the value if necessary.
    if (.not. associated(p_rvalue%p_Iarray)) then
      allocate(p_rvalue%p_Iarray(size(value)))
    else if (size(p_rvalue%p_Iarray) .ne. size(value)) then
      deallocate(p_rvalue%p_Iarray)
      allocate(p_rvalue%p_Iarray(size(value)))
    end if

    ! Set the value
    p_rvalue%p_Iarray = value

  end subroutine collct_setvalue_intarr

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_iarrp (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
  !
  ! The routine creates a new memory block on the heap and copies the
  ! data to there. There is no reference to the old array stored.
  !
  ! The routine saves a pointer to an integer array.
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  integer, dimension(:), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGERARRP,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    p_rvalue%p_Iarray => value

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_real (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  real(DP), intent(in) :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_REAL,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%dvalue = value

  end subroutine collct_setvalue_real

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_realarr (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
  !
  ! The routine creates a new memory block on the heap and copies the
  ! data to there. There is no reference to the old array stored.
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  real(DP), dimension(:),intent(in) :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_REALARR,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! (Re-)allocate memory for the value if necessary.
    if (.not. associated(p_rvalue%p_Darray)) then
      allocate(p_rvalue%p_Darray(size(value)))
    else if (size(p_rvalue%p_Darray) .ne. size(value)) then
      deallocate(p_rvalue%p_Darray)
      allocate(p_rvalue%p_Darray(size(value)))
    end if

    ! Set the value
    p_rvalue%p_Darray = value

  end subroutine collct_setvalue_realarr

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_rarrp (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Sets the value of the parameter: sparameter=value
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
  !
  ! The routine creates a new memory block on the heap and copies the
  ! data to there. There is no reference to the old array stored.
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  real(DP), dimension(:),intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_REALARRP,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    p_rvalue%p_Darray => value

  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_discr (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_spatialDiscretisation), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCR,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rdiscretisation => value

  end subroutine collct_setvalue_discr

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_bldiscr (rcollection, sparameter, value, badd, &
                                      ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_blockDiscretisation), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BLDISCR,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rblockDiscretisation => value

  end subroutine collct_setvalue_bldiscr

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_tria (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_triangulation), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_TRIA,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rtriangulation => value

  end subroutine collct_setvalue_tria

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_bdry (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_boundary), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARY,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rboundary => value

  end subroutine collct_setvalue_bdry

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_bdreg (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_boundaryRegion), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BDREGION,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rboundaryRegion => value

  end subroutine collct_setvalue_bdreg

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_bc (rcollection, sparameter, value, badd, &
                                 ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_boundaryConditions), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARYCOND,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rboundaryConditions => value

  end subroutine collct_setvalue_bc

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_vecsca (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_vectorScalar), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAVECTOR,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rvectorScalar => value

  end subroutine collct_setvalue_vecsca

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_matsca (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_matrixScalar), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAMATRIX,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rmatrixScalar => value

  end subroutine collct_setvalue_matsca

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_vec (rcollection, sparameter, value, badd, &
                                  ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_vectorBlock), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKVECTOR,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rvector => value

  end subroutine collct_setvalue_vec

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_mat (rcollection, sparameter, value, badd, &
                                  ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_matrixBlock), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKMATRIX,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rmatrix => value

  end subroutine collct_setvalue_mat

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_parlst (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_parlist), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_PARAMETERS,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rparlist => value

  end subroutine collct_setvalue_parlst

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_linsol (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_linsolNode), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_LINSOL,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rlinearSolver => value

  end subroutine collct_setvalue_linsol

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_discbc (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_discreteBC), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCRBC,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rdiscreteBC => value

  end subroutine collct_setvalue_discbc

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_ilvpsc (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_interlevelProjectionScalar), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJSC,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rilvProjectionSc => value

  end subroutine collct_setvalue_ilvpsc

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_ilvp (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_interlevelProjectionBlock), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJ,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rilvProjection => value

  end subroutine collct_setvalue_ilvp

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_coll (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_collection), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_COLLECTION,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rcollection => value

  end subroutine collct_setvalue_coll

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_pars (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_fparser), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FPARSER,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rparser => value

  end subroutine collct_setvalue_pars

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_geom (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_geometryObject), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GEOMETRY,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rgeometry => value

  end subroutine collct_setvalue_geom

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_fchn (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_filterChain), dimension(:), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FILTERCHAIN,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_RfilterChain => value

  end subroutine collct_setvalue_fchn

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_hadapt (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_hadapt), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_HADAPT,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rhadapt => value

  end subroutine collct_setvalue_hadapt

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_afcstab (rcollection, sparameter, value, badd, &
                                      ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_afcstab), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_AFCSTAB,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rafcstab => value

  end subroutine collct_setvalue_afcstab

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_timer (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_timer), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_TIMER,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rtimer => value

  end subroutine collct_setvalue_timer

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_mshh (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_meshHierarchy), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_MSHHIERARCHY,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rmeshHierarchy => value

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_fesp (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_feSpaceLevel), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FESPACE,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rfeSpace => value

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_feh (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_feHierarchy), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_FEHIERARCHY,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rfeHierarchy => value

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_tsh (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_timescaleHierarchy), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_TSHIERARCHY,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rtsHierarchy => value

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_mlprjh (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_interlevelProjectionHier), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_MLPRJHIERARCHY,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rmlprjHierarchy => value

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_graph (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_graph), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GRAPH,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rgraph => value

  end subroutine collct_setvalue_graph

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_particles(rcollection, sparameter, value, badd, &
                                       ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_particleCollection), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_PARTICLES,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rparticles => value

  end subroutine collct_setvalue_particles

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_particles3D(rcollection, sparameter, value, badd, &
                                         ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_particleCollection3D), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_PARTICLES3D,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rparticles3D => value

  end subroutine collct_setvalue_particles3D

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_gfemset(rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_groupFEMSet), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GROUPFEMSET,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rgroupFEMSet => value

  end subroutine collct_setvalue_gfemset

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_gfemblk(rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_groupFEMBlock), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GROUPFEMBLOCK,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rgroupFEMBlock => value

  end subroutine collct_setvalue_gfemblk

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_cubinfo(rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_scalarCubatureInfo), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_CUBINFO,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rcubatureInfo => value

  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine collct_setvalue_genobj(rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName)
!<description>
  ! Stores a pointer to 'value' using the parameter name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>

!<inputoutput>

  ! The parameter list.
  type(t_collection), intent(inout) :: rcollection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter.
  type(t_genericObject), intent(in), target :: value

  ! Whether to add the variable if it does not exist.
  ! =false: do not add the variable, throw an error
  ! =true : add the variable
  logical, intent(in) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  integer, intent(in), optional :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  character(LEN=*), intent(in), optional :: ssectionName

!</input>

!</subroutine>

    ! local variables
    type(t_collctValue), pointer :: p_rvalue
    logical :: bexists

    ! Get the pointer to the parameter. Add the parameter if necessary
    call collct_getvalue_struc (rcollection, sparameter, COLLCT_GENOBJ,&
                                badd, p_rvalue, ilevel, bexists, ssectionName)

    ! Set the value
    p_rvalue%p_rgenericObject => value

  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine collct_copyQuickAccess(rcollectionSrc, rcollectionDest)

!<description>
  ! This subroutine copies all quick access arrays from the collection
  ! structure rcollectionSrc to the collection structure rcollectionDest
!</description>

!<input>
  ! Source collection structure
  type(t_collection), intent(in) :: rcollectionSrc
!</input>

!<inputoutput>
  ! Destination collection structure
  type(t_collection), intent(inout) :: rcollectionDest
!</inputoutput>

!</subroutine>

    rcollectionDest%IquickAccess = rcollectionSrc%IquickAccess
    rcollectionDest%DquickAccess = rcollectionSrc%DquickAccess
    rcollectionDest%SquickAccess = rcollectionSrc%SquickAccess

    rcollectionDest%p_rvectorQuickAccess1 => rcollectionSrc%p_rvectorQuickAccess1
    rcollectionDest%p_rvectorQuickAccess2 => rcollectionSrc%p_rvectorQuickAccess2
    rcollectionDest%p_rvectorQuickAccess3 => rcollectionSrc%p_rvectorQuickAccess3
    rcollectionDest%p_rvectorQuickAccess4 => rcollectionSrc%p_rvectorQuickAccess4

    rcollectionDest%p_rmatrixQuickAccess1 => rcollectionSrc%p_rmatrixQuickAccess1
    rcollectionDest%p_rmatrixQuickAccess2 => rcollectionSrc%p_rmatrixQuickAccess2
    rcollectionDest%p_rmatrixQuickAccess3 => rcollectionSrc%p_rmatrixQuickAccess3
    rcollectionDest%p_rmatrixQuickAccess4 => rcollectionSrc%p_rmatrixQuickAccess4

    rcollectionDest%p_rparlistQuickAccess1 => rcollectionSrc%p_rparlistQuickAccess1
    rcollectionDest%p_rparlistQuickAccess2 => rcollectionSrc%p_rparlistQuickAccess2
    rcollectionDest%p_rparlistQuickAccess3 => rcollectionSrc%p_rparlistQuickAccess3
    rcollectionDest%p_rparlistQuickAccess4 => rcollectionSrc%p_rparlistQuickAccess4

  end subroutine collct_copyQuickAccess

end module collection
