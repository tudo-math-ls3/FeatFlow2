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
!#
!# The list supports 'level tags' and 'section tags" to group values:
!# 
!# 1.) There's one unnamed section ('') in the collection for general purpose
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
!# Additionally to these 'named variables', the collection contains
!# a small integer and double precision 'temporary quick access'
!# array directly in the t_collection structure with size of at
!# least 128 elemens. 
!# Care must be taken when an application wants to use these arays.
!# Not maintained by the collection internally, they are of 
!# *pure temprary nature*. Long-life information should never be stored
!# in there! They can be used e.g. in assembly routines to provide 
!# callback routines with information that must quickly be accessed 
!# without the need of asking the collection for a named variable.
!# Their content can be described at best as 'only valid for one
!# operation' (e.g. for one matrix asembly).
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
!# GUILDLINE
!# ---------
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
!#  The collection is a black-box container. It's designed to hold low- and 
!#  mid-level information (numbers, strings, as well as linear algebra stuff
!#  like vectors, matrices) that is used for solving linear problems. It's
!#  clearly NOT designed that a user adds problem-specific data structures
!#  there. Trying this will result in circular dependencies! Adding information
!#  from a problem-dependent structure to the collection with the 'add' 
!#  routines is completely right, but don't extent the collection
!#  by new data types / structures unless you are absolutely sure that this 
!#  is relevant for all users and that you won't get circular dependencies!
!#
!# USE THE COLLECTION IS AN INDEX
!#
!#  Although it's possible to use the collection as a data container for most
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
!#  but don't use this feature for storing long-life information!
!#  Example:
!#
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
!#  It's possible to store integer / real arrays in the collection. These
!#  should only be small and the application must keep track of the size
!#  of these arrays. To store larger arrays, use the STORAGE module to
!#  allocate memory and store the handle of the memory block in the collection!
!# 
!# </purpose>
!##############################################################################

MODULE collection

  USE fsystem
  USE io
  USE linearsystemscalar
  USE linearsystemblock
  USE paramlist
  USE spatialdiscretisation
  USE linearsolver
  USE boundary
  USE multilevelprojection
  USE filtersupport
  USE fparser
  USE genoutput
  USE geometry
  USE hadaptaux
  USE afcstabilisation
  
  IMPLICIT NONE

!<constants>

!<constantblock>
  
  ! Maximum length of a section name.
  INTEGER, PARAMETER :: COLLCT_MLSECTION = 64

  ! Maximum length of the name of a value: 32 characters
  INTEGER, PARAMETER :: COLLCT_MLNAME = 32

  ! Minimum number of free 'slots' per 'level' subgroup.
  ! If there are too many values in a 'level' subgroup, the
  ! structure is dynamically extended in terms of collct_NVALUES
  ! entries.
  INTEGER, PARAMETER :: COLLCT_NVALUES = 32

  ! Minimum number of level subgroups (additionally to the 0th).
  ! If there are too many levels in a parameter block, the
  ! structure is dynamically extended in terms of COLLCT_NLEVELS
  ! entries.
  INTEGER, PARAMETER :: COLLCT_NLEVELS = 9

  ! Minimum number of sections - additionally to the unnamed one.
  ! If there are too many sections in a collection block, the
  ! structure is dynamically extended in terms of COLLCT_NSECTIONS
  ! entries.
  INTEGER, PARAMETER :: COLLCT_NSECTIONS = 5

  ! Length of each 'quick-access' array in the collection.
  INTEGER, PARAMETER :: COLLCT_QALENGTH = 128
  
  ! Number of 'quick-access' strings in the collection
  INTEGER, PARAMETER :: COLLCT_QASTRINGS = 4
!</constantblock>

!<constantblock description="Type identifier for values">
  
  ! Undefined value
  INTEGER, PARAMETER :: COLLCT_UNDEFINED = 0
  
  ! Character value
  INTEGER, PARAMETER :: COLLCT_CHARACTER = 1

  ! String value
  INTEGER, PARAMETER :: COLLCT_STRING    = 2
  
  ! Integer value
  INTEGER, PARAMETER :: COLLCT_INTEGER   = 3
  
  ! Double precision REAL value 
  INTEGER, PARAMETER :: COLLCT_REAL      = 4

  ! Scalar discretisation structure
  INTEGER, PARAMETER :: COLLCT_DISCR     = 5

  ! Block discretisation structure
  INTEGER, PARAMETER :: COLLCT_BLDISCR   = 5

  ! Triangulation structure
  INTEGER, PARAMETER :: COLLCT_TRIA      = 7

  ! A scalar vector
  INTEGER, PARAMETER :: COLLCT_SCAVECTOR = 8

  ! A scalar matrix
  INTEGER, PARAMETER :: COLLCT_SCAMATRIX = 9

  ! A block vector
  INTEGER, PARAMETER :: COLLCT_BLKVECTOR = 10

  ! A block matrix
  INTEGER, PARAMETER :: COLLCT_BLKMATRIX = 11
  
  ! A parameter structure
  INTEGER, PARAMETER :: COLLCT_PARAMETERS = 12

  ! A linear solver structure
  INTEGER, PARAMETER :: COLLCT_LINSOL     = 13

  ! Discretised-boundary-conditions structure
  INTEGER, PARAMETER :: COLLCT_DISCRBC    = 14

  ! A domain
  INTEGER, PARAMETER :: COLLCT_BOUNDARY   = 15

  ! Scalar analytic boundary conditions
  INTEGER, PARAMETER :: COLLCT_BOUNDARYCOND = 16

  ! Scalar interlevel projection structure
  INTEGER, PARAMETER :: COLLCT_INTERLVPRJSC = 17
  
  ! Block interlevel projection structure
  INTEGER, PARAMETER :: COLLCT_INTERLVPRJ   = 18

  ! A parser structure for a 'precompiled' expression
  INTEGER, PARAMETER :: COLLCT_FPARSER      = 19

  ! A filter chain
  INTEGER, PARAMETER :: COLLCT_FILTERCHAIN  = 20

  ! Integer value type array
  INTEGER, PARAMETER :: COLLCT_INTEGERARR   = 21
  
  ! Double precision REAL value type array
  INTEGER, PARAMETER :: COLLCT_REALARR      = 22

  ! Geometry object
  INTEGER, PARAMETER :: COLLCT_GEOMETRY     = 24

  ! The collection structure itself
  INTEGER, PARAMETER :: COLLCT_COLLECTION   = 25

  ! Grid adaptation object
  INTEGER, PARAMETER :: COLLCT_HADAPT       = 26

  ! Stabilisation structure
  INTEGER, PARAMETER :: COLLCT_AFCSTAB      = 27

!</constantblock>

!</constants>

!<types>
  
!<typeblock>
  
  ! This structure realises a value in the scalar collection.
  ! It contains a type and some data variables that contain the value.

  TYPE t_collctValue
  
    PRIVATE
   
    ! The type of the section (COLLCT_UNDEFINED, COLLCT_CHARACTER,...)
    INTEGER :: itype = COLLCT_UNDEFINED
    
    ! Name of the value
    CHARACTER(LEN=COLLCT_MLNAME) :: sname
    
    ! Character value
    CHARACTER :: csvalue = " "

    ! String value
    CHARACTER, DIMENSION(:), POINTER :: p_svalue => NULL()
    
    ! Integer value
    INTEGER(I32) :: ivalue = 0
    
    ! Double precision value 
    REAL(DP)     :: dvalue = 0.0_DP

    ! Pointer to a spatial discretisation structure
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation => NULL()

    ! Pointer to a block discretisation structure
    TYPE(t_blockDiscretisation), POINTER :: p_rblockDiscretisation => NULL()

    ! Pointer to a triangulation structure, 2D
    TYPE(t_triangulation), POINTER :: p_rtriangulation => NULL()

    ! Pointer to a scalar vector
    TYPE(t_vectorScalar), POINTER :: p_rvectorScalar => NULL()

    ! Pointer to a scalar matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixScalar => NULL()
    
    ! Pointer to a block vector
    TYPE(t_vectorBlock), POINTER :: p_rvector => NULL()

    ! Pointer to a block matrix
    TYPE(t_matrixBlock), POINTER :: p_rmatrix => NULL()

    ! Pointer to a parameter list
    TYPE(t_parlist), POINTER     :: p_rparlist => NULL()

    ! Pointer to a linear solver structure
    TYPE(t_linsolNode), POINTER  :: p_rlinearSolver => NULL()

    ! A structure containing discretised boundary conditions
    TYPE(t_discreteBC), POINTER  :: p_rdiscreteBC     => NULL()

    ! Pointer to a domain
    TYPE(t_boundary), POINTER      :: p_rboundary => NULL()

    ! Pointer to a geometry object
    TYPE(t_geometryObject), POINTER :: p_rgeometry => NULL()

    ! Pointer to scalar boundary conditions
    TYPE(t_boundaryConditions), POINTER      :: p_rboundaryConditions => NULL()

    ! Pointer to a scalar interlevel projection structure
    TYPE(t_interlevelProjectionScalar), POINTER :: p_rilvprojectionSc => NULL()

    ! Pointer to a scalar interlevel projection structure
    TYPE(t_interlevelProjectionBlock), POINTER  :: p_rilvprojection => NULL()

    ! Pointer to a parser object
    TYPE(t_fparser), POINTER                    :: p_rparser => NULL()

    ! Pointer to a filter chain
    TYPE(t_filterChain), DIMENSION(:), POINTER  :: p_RfilterChain => NULL()

    ! Pointer to an integer precision array
    INTEGER(I32), DIMENSION(:), POINTER         :: p_Iarray => NULL()

    ! Pointer to an real precision array
    REAL(DP), DIMENSION(:), POINTER             :: p_Darray => NULL()

    ! Pointer to a collection structure
    TYPE(t_collection), POINTER                 :: p_rcollection => NULL()

    ! Pointer to a grid adaptation structure
    TYPE(t_hadapt), POINTER                     :: p_rhadapt => NULL()

    ! Pointer to a stabilisation structure
    TYPE(t_afcstab), POINTER                    :: p_rafcstab => NULL()
  END TYPE
  
!</typeblock>
  
!<typeblock>
  
  ! This structure realises a level structure containing multiple values.
  
  TYPE t_collctLevel
  
    PRIVATE
   
    ! Actual number of values in this section.
    INTEGER :: ivalueCount = 0
    
    ! Indicates whether the first ivalueCount values are 'en-block'
    ! or if there are 'holes' inbetween because of deleted values.
    LOGICAL :: bisFull = .TRUE.
    
    ! A list of values.
    TYPE(t_collctValue), DIMENSION(:), POINTER :: p_Rvalues => NULL()
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! This structure realises a section. Each section contains a list of
  ! level substructures and one substructure not assigned to a level,
  ! containing level independent parameters.
  
  TYPE t_collctSection
  
    PRIVATE
  
    ! The name of the section. '' identifies an/the unnamed section.
    CHARACTER(LEN=COLLCT_MLSECTION) :: ssectionName = ''

    ! Actual number of levels in this section. There's at least
    ! one section - the unnamed section. If this value is =0, the parameter
    ! list is currently not initialised.
    INTEGER :: ilevelCount = 0
    
    ! The 0th level collecting general variables, level-independent
    TYPE(t_collctLevel) :: rlevel0
    
    ! A list of levels or NULL(), if there are no additional 'level'
    ! subgroups than level 0.
    TYPE(t_collctLevel), DIMENSION(:), POINTER :: p_Rlevels => NULL()
    
  END TYPE
  
!<typeblock>
  
!</typeblock>
  
  ! This structure realises a collection.
  ! A collection contains one unnamed section (containing general data)
  ! and arbitrary many named sections.
  ! Each section contains a level structure for level 0 (containing
  ! level independent data) and a number of levels where data can
  ! be stored.
  
  TYPE t_collection
  
    ! USER DEFINED:
    ! Quick-access array for integers. This is a short, temporary
    ! array that can be used by the application to save intermediate
    ! values (e.g. some information that is to pass and to be accessed
    ! by callback routines). 
    ! No long life information should be stored here. Long-life
    ! information should be named and added as such to the collection 
    ! directly. The collection itself does not maintain this array!
    INTEGER(I32), DIMENSION(COLLCT_QALENGTH) :: IquickAccess = 0
  
    ! USER DEFINED:
    ! Quick-access array for reals. This is a short, temporary
    ! array that can be used by the application to save intermediate
    ! values (e.g. some information that is to pass and to be accessed
    ! by callback routines). 
    ! No long life information should be stored here. Long-life
    ! information should be named and added as such to the collection 
    ! directly. The collection itself does not maintain this array!
    REAL(DP), DIMENSION(COLLCT_QALENGTH) :: DquickAccess = 0.0_DP
  
    ! USER DEFINED:
    ! A simple string array for direct access.
    CHARACTER(LEN=SYS_STRLEN), DIMENSION(COLLCT_QASTRINGS) :: SquickAccess

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    TYPE(t_vectorBlock), POINTER :: p_rvectorQuickAccess1 => NULL()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    TYPE(t_vectorBlock), POINTER :: p_rvectorQuickAccess2 => NULL()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    TYPE(t_vectorBlock), POINTER :: p_rvectorQuickAccess3 => NULL()

    ! USER DEFINED:
    ! A quick access pointer to a vector.
    TYPE(t_vectorBlock), POINTER :: p_rvectorQuickAccess4 => NULL()
    
    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    TYPE(t_vectorBlock), POINTER :: p_rmatrixQuickAccess1 => NULL()

    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    TYPE(t_vectorBlock), POINTER :: p_rmatrixQuickAccess2 => NULL()

    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    TYPE(t_vectorBlock), POINTER :: p_rmatrixQuickAccess3 => NULL()

    ! USER DEFINED:
    ! A quick access pointer to a matrix.
    TYPE(t_vectorBlock), POINTER :: p_rmatrixQuickAccess4 => NULL()

    ! USER DEFINED:
    ! A quick access pointer to another collection structure.
    ! Can be used to link collection structures.
    TYPE(t_collection), POINTER :: p_rnextCollection => NULL()
    
    ! Actual number of sections in this collection.
    ! This is at least one, as every collection contains at least an
    ! unnamed section.
    INTEGER :: isectionCount = 0
    
    ! A list of levels or NULL(), if there are no additional 'level'
    ! subgroups than level 0.
    TYPE(t_collctSection), DIMENSION(:), POINTER :: p_Rsections => NULL()
    
  END TYPE

!</typeblock>

!</types>

  PRIVATE :: collct_initlevel, collct_donelevel, collct_realloclevel
  PRIVATE :: collct_initsection, collct_donesection, collct_reallocsection
  PRIVATE :: collct_realloccollection
  PRIVATE :: collct_fetchparameter_indir, collct_fetchparameter_direct
  
  PRIVATE :: collct_addlevel_indir, collct_queryvalue_indir, collct_addvalue
  PRIVATE :: collct_fetchsection, collct_fetchlevel,collct_getmaxlevel_indir 
  
  INTERFACE collct_addlevel
    MODULE PROCEDURE collct_addlevel_indir
    MODULE PROCEDURE collct_addlevel_direct
    MODULE PROCEDURE collct_addlevel_all
  END INTERFACE
  
  INTERFACE collct_queryvalue
    MODULE PROCEDURE collct_queryvalue_direct 
    MODULE PROCEDURE collct_queryvalue_indir
  END INTERFACE
  
CONTAINS
  
  ! ***************************************************************************

  ! Internal subroutine: Initialise a newly created level subgroup.
  
  SUBROUTINE collct_initlevel (rcollctLevel)
  
  TYPE(t_collctLevel), INTENT(INOUT) :: rcollctLevel
  
  ! Nullify the level pointer. Is allocated when the first entry is put
  ! into the list.
  NULLIFY(rcollctLevel%p_Rvalues)
  
  ! No variables in here for the moment
  rcollctLevel%ivalueCount = 0
  
  ! No holes.
  rcollctLevel%bisFull = .TRUE.
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Releases a level and all values inside from memory.
  
  SUBROUTINE collct_donelevel (rcollctLevel)
  
  TYPE(t_collctLevel), INTENT(INOUT) :: rcollctLevel
  
  ! Deallocate the value list on the current level if there is one.
  IF (ASSOCIATED(rcollctLevel%p_Rvalues)) THEN
    DEALLOCATE(rcollctLevel%p_Rvalues)
  END IF
  
  ! No variables in here for the moment
  rcollctLevel%ivalueCount = 0
  
  ! No holes anymore.
  rcollctLevel%bisFull = .TRUE.
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Reallocate a level.
  ! This increases the size of a level subgroup by reallocation of the
  ! arrays. 
  
  SUBROUTINE collct_realloclevel (rcollctLevel, inewsize)
  
  ! The section to reallocate.
  TYPE(t_collctLevel), INTENT(INOUT) :: rcollctLevel
  
  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the level should be able to handle.
  INTEGER, INTENT(IN) :: inewsize
  
  ! local variables
  INTEGER :: sz

  ! Pointers to new lists for replacing the old.
  TYPE(t_collctValue), DIMENSION(:), POINTER :: p_Rvalues
  
  sz = MAX(SIZE(rcollctLevel%p_Rvalues),inewsize)

  IF (SIZE(rcollctLevel%p_Rvalues) .EQ. sz) RETURN ! nothing to do
  
  ! Allocate the pointers for the new list
  ALLOCATE(p_Rvalues(sz))
  
  ! Copy the content of the old ones
  p_Rvalues(1:sz) = rcollctLevel%p_Rvalues (1:sz)
  
  ! Throw away the old array, replace by the new one
  DEALLOCATE(rcollctLevel%p_Rvalues)
  
  rcollctLevel%p_Rvalues => p_Rvalues
  
  END SUBROUTINE

  
  ! ***************************************************************************

  ! Internal subroutine: Add a new value structure to a level and return
  ! a pointer to it. Returns NULL() if the section or the level does not
  ! exist.
  
  SUBROUTINE collct_addvalue (rcollection, ssectionName, sparamName, itype, &
                              ilevel, p_rvalue)
  
  ! The collection
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The name of the section; '' identifies the unnamed section
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The name of the value to add. Must be <> ''!
  CHARACTER(LEN=*), INTENT(IN) :: sparamName
  
  ! The type of the parameter
  INTEGER, INTENT(IN) :: itype

  ! The level where to add; 0=level independent part
  INTEGER :: ilevel
  
  ! Output: The pointer of the newly created value
  TYPE(t_collctValue), POINTER :: p_rvalue
  
  ! local variables
  TYPE(t_collctSection), POINTER :: p_rSection
  TYPE(t_collctLevel), POINTER :: p_rLevel
  TYPE(t_collctValue), DIMENSION(:), POINTER :: p_Rvalues
  INTEGER :: i
  
  NULLIFY(p_rvalue)
  
  ! The name must be given!
  IF (sparamName .EQ. '') RETURN
  
  ! Get the section
  CALL collct_fetchsection(rcollection, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) RETURN
  
  ! Get the level
  i=MAX(ilevel,0)
  IF (i .GT. p_rsection%ilevelCount) RETURN
  
  IF (i .EQ. 0) THEN
    p_rLevel => p_rsection%rlevel0
  ELSE
    p_rLevel => p_rsection%p_Rlevels(i)
  END IF
  
  ! Get space for the new structure in the list
  
  IF (.NOT. ASSOCIATED(p_rLevel%p_Rvalues)) THEN
  
    ! Create an array for the entries
    ALLOCATE(p_rLevel%p_Rvalues(COLLCT_NVALUES))
    
    ! Take a pointer to the first position for storing data.
    p_rvalue => p_Rlevel%p_Rvalues(1)
    
  ELSE IF (p_rLevel%ivalueCount .GE. SIZE(p_rLevel%p_Rvalues)) THEN
  
    ! Reallocate the entry array to get space for the new entry
    ALLOCATE(p_Rvalues(SIZE(p_rLevel%p_Rvalues)+COLLCT_NVALUES))
    p_Rvalues(1:p_rLevel%ivalueCount) = p_rLevel%p_Rvalues (1:p_rLevel%ivalueCount)
    DEALLOCATE(p_rLevel%p_Rvalues)
    p_rLevel%p_Rvalues => p_Rvalues
    
    ! Store the value at the new free position.
    p_rvalue => p_Rlevel%p_Rvalues(p_rLevel%ivalueCount+1)
    
    ! There are definitely no holes anymore in the list on this
    ! level, as all free positions are filled up when we are here.
    p_rLevel%bisFull = .TRUE.
    
  ELSE
    ! Are there holes or can we directly take the maximum position?
    IF (p_rLevel%bisFull) THEN
      ! No holes, take the next free position.
      i = p_rLevel%ivalueCount + 1
    ELSE
      ! Find an empty position in the array.
      DO i=1,SIZE(p_rLevel%p_Rvalues)
        IF (p_rLevel%p_Rvalues(i)%itype .EQ. COLLCT_UNDEFINED) EXIT
      END DO
    
      ! This must work, otherwise ivalueCount is wrong!
      IF (i .GT. SIZE(p_rLevel%p_Rvalues)) THEN
        PRINT *,'Error in collct_addvalue: Collection structure inconsistent!'
        CALL sys_halt()
      END IF
      
      ! If i is larger than ivalueCount, all holes are filled up again,
      ! i.e. we add 'behind' all elements.
      IF (i .GT. p_rLevel%ivalueCount) p_rLevel%bisFull = .TRUE.
      
    END IF

    ! Take that value for storing data.
    p_rvalue => p_Rlevel%p_Rvalues(i)
    
  END IF

  ! Fill the value with initial data and increase the counter.
  p_rLevel%ivalueCount = p_rLevel%ivalueCount + 1
  p_rvalue%itype = itype

  ! Store the name in upper case
  p_rvalue%sname = sys_upcase(sparamName)
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Initialise a newly created section
  
  SUBROUTINE collct_initsection (rcollctSection, ssectionName)
  
  ! The collection
  TYPE(t_collctSection), INTENT(INOUT) :: rcollctSection
  
  ! The name for the new section
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName
  
  ! Simply allocate the pointers with an empty list
  ALLOCATE(rcollctSection%p_Rlevels(COLLCT_NVALUES))
  
  ! No variables in here for the moment
  rcollctSection%ilevelCount = 0
  
  ! Set the section name - uppercase
  rcollctSection%ssectionName = ssectionName
  CALL sys_toupper(rcollctSection%ssectionName)
  
  ! Init level 0
  CALL collct_initlevel (rcollctSection%rlevel0)
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Releases a section structure from memory.
  
  SUBROUTINE collct_donesection (rcollctSection)
  
  TYPE(t_collctSection), INTENT(INOUT) :: rcollctSection
  
  INTEGER :: i
  
  ! Clean up level 0
  CALL collct_donelevel (rcollctSection%rlevel0)
  
  ! Clean up the other levels in this section
  DO i=rcollctSection%ilevelCount,1,-1
    CALL collct_donelevel (rcollctSection%p_Rlevels(i))
  END DO
  
  ! Deallocate the pointers with an empty list
  DEALLOCATE(rcollctSection%p_Rlevels)
  
  ! No variables in here for the moment
  rcollctSection%ilevelCount = 0
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Reallocate a section.
  ! This increases the size of a section (i.e. the number of levels that
  ! can be stored in a section) by reallocation of the arrays.
  
  SUBROUTINE collct_reallocsection (rcollctSection, inewsize)
  
  ! The section to reallocate.
  TYPE(t_collctSection), INTENT(INOUT) :: rcollctSection
  
  ! The new 'size' of the section, i.e. the new number of levels,
  ! the section should be able to handle.
  INTEGER, INTENT(IN) :: inewsize
  
  ! local variables
  INTEGER :: sz

  ! Pointers to new lists for replacing the old.
  TYPE(t_collctLevel), DIMENSION(:), POINTER :: p_Rlevels
  
  sz = MAX(SIZE(rcollctSection%p_Rlevels),inewsize)

  IF (SIZE(rcollctSection%p_Rlevels) .EQ. sz) RETURN ! nothing to do
  
  ! Allocate the pointers for the new list
  ALLOCATE(p_Rlevels(sz))
  
  ! Copy the content of the old ones
  p_Rlevels(1:sz) = rcollctSection%p_Rlevels (1:sz)
  
  ! Throw away the old array, replace by the new one
  DEALLOCATE(rcollctSection%p_Rlevels)
  
  rcollctSection%p_Rlevels => p_Rlevels
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Reallocate the collection.
  ! This increases the size of a a collection (i.e. the number of sections that
  ! can be stored in a collection) by reallocation of the arrays.
  
  SUBROUTINE collct_realloccollection (rcollection, inewsize)
  
  ! The section list to reallocate.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the section should be able to handle.
  INTEGER, INTENT(IN) :: inewsize
  
  ! local variables
  
  INTEGER :: sz

  ! Pointers to new lists for replacing the old.
  TYPE(t_collctSection), DIMENSION(:), POINTER :: p_Rsections
  
  sz = MAX(SIZE(rcollection%p_Rsections),inewsize)

  IF (SIZE(rcollection%p_Rsections) .EQ. sz) RETURN ! nothing to do

  ! Allocate the pointers for the new lists
  ALLOCATE(p_Rsections(sz))
  
  ! Copy the content of the old ones
  p_Rsections(1:sz) = rcollection%p_Rsections (1:sz)
  
  ! Throw away the old array, replace by the new one
  DEALLOCATE(rcollection%p_Rsections)
  
  rcollection%p_Rsections => p_Rsections
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
  ! Internal subroutine: Search in a level subgroup for a parameter
  ! and return the index - or 0 if the parameter does not exist.

  SUBROUTINE collct_fetchparameter_indir (rlevel, sname, iparamnum) 

  ! The level where to search
  TYPE(t_collctLevel), INTENT(IN) :: rlevel
  
  ! The parameter name to look for. 
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! The number of the parameter in the list or 0 if it does not exist.
  INTEGER, INTENT(OUT) :: iparamnum
  
  ! local variables
  INTEGER :: i,nsections
  CHARACTER(LEN=COLLCT_MLNAME)  :: sname2
  
  ! Convert the name to uppercase.
  sname2 = sname
  CALL sys_toupper (sname2)
  
  iparamnum = 0
  
  ! If the parameter list is empty, the section does not exist for sure
  IF (rlevel%ivalueCount .EQ. 0) RETURN
  
  ! Loop through all sections to see if the section exists
  IF (rlevel%bisfull) THEN
    nsections = rlevel%ivalueCount
  ELSE
    nsections = SIZE(rlevel%p_Rvalues)
  END IF

  DO i=1,nsections
    IF (rlevel%p_Rvalues(i)%sname .EQ. sname2) THEN
      iparamnum = i
      RETURN
    END IF
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
  ! Internal subroutine: Search in a level subgroup for a parameter
  ! and return a pointer to it - or NULL if the parameter does not exist

  SUBROUTINE collct_fetchparameter_direct (rcollection, ssectionName, ilevel, &
             sparameter, p_rvalue) 

  ! The collection where to search
  TYPE(t_collection), INTENT(IN) :: rcollection
  
  ! The section name where to search - '' identifies the unnamed section
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName
  
  ! The level where to search; maybe 0 for level-independent parameters
  INTEGER, INTENT(IN) :: ilevel
  
  ! The parameter name to look for. 
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! A pointer to the parameter.
  TYPE(t_collctValue), POINTER :: p_rvalue
  
  ! local variables
  INTEGER :: i,ilv
  TYPE(t_collctSection), POINTER :: p_rsection
  TYPE(t_collctLevel), POINTER :: p_rlevel
  
  ! Some basic checks
  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    CALL sys_halt()
  END IF

  NULLIFY(p_rvalue)

  CALL collct_fetchsection(rcollection, ssectionName, p_rsection) 
  
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    ! Section does not exist - return NULL
    RETURN
  END IF
  
  ! Get the level - if it exists.
  ilv = MAX(0,ilevel)
  IF(ilv .GT. p_rsection%ilevelCount) THEN
    ! Level does not exist - return NULL
    RETURN
  END IF
  
  IF (ilv .EQ. 0) THEN
    p_rlevel => p_rsection%rlevel0
  ELSE
    p_rlevel => p_rsection%p_Rlevels(ilv)
  END IF
  
  ! Get the index of the value
  CALL collct_fetchparameter_indir (p_rlevel, sparameter, i) 
  
  ! and finally the pointer to it
  IF (i .NE. 0) THEN
    p_rvalue => p_rlevel%p_Rvalues(i)
  END IF

  END SUBROUTINE

  ! ***************************************************************************
  
  ! Internal subroutine: Search in a collection for a section
  ! and return a pointer to the section - or NULL() if the section does 
  ! not exist.

  SUBROUTINE collct_fetchsection(rcollection, sname, p_rsection, isectionIndex) 

  ! The section.
  TYPE(t_collection), INTENT(IN) :: rcollection
  
  ! The parameter name to look for. 
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! A pointer to the section or NULL() if it does not exist.
  TYPE(t_collctSection), POINTER :: p_rsection
  
  ! Optional output: Index of the section in the section array.
  INTEGER, INTENT(OUT), OPTIONAL :: isectionIndex
  
  ! local variables
  INTEGER :: i
  CHARACTER(LEN=COLLCT_MLSECTION) :: sname2
  
  ! Convert the name to uppercase.
  sname2 = sname
  CALL sys_toupper (sname2)
  
  NULLIFY(p_rsection)
  
  ! If the parameter list is empty, the section does not exist for sure
  IF (rcollection%isectionCount .EQ. 0) RETURN
  
  ! Loop through all sections to see if the section exists
  DO i=1,rcollection%isectionCount
    IF (rcollection%p_Rsections(i)%ssectionName .EQ. sname2) THEN
      p_rsection => rcollection%p_Rsections(i)
      
      IF (PRESENT(isectionIndex)) isectionIndex = i
      
      RETURN
    END IF
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
  ! Internal subroutine: Search in a collection for a level in a section
  ! and return a pointer to the level - or NULL)( if the level or the section
  ! does not exist).

  SUBROUTINE collct_fetchlevel(rcollection, ssectionName, ilevel, p_rlevel) 

  ! The section.
  TYPE(t_collection), INTENT(IN) :: rcollection
  
  ! The section name to look for. 
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName
  
  ! The level number where to search for
  INTEGER, INTENT(IN) :: ilevel
  
  ! OUTPUT: A pointer to the level - or NULL() if the level does not exist
  TYPE(t_collctLevel), POINTER :: p_rlevel
  
  ! local variables
  INTEGER :: i
  TYPE(t_collctSection), POINTER :: p_rsection
  
  NULLIFY(p_rlevel)
  
  ! Get the section
  CALL collct_fetchsection(rcollection, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) RETURN
  
  ! Get the level
  i = MAX(0,ilevel)
  IF (i .GT. p_rsection%ilevelCount) RETURN
  
  IF (i .EQ. 0) THEN
    ! Get level 0 - the always-existing level
    p_rlevel => p_rsection%rlevel0
  ELSE
    p_rlevel => p_rsection%p_Rlevels(i)
  END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE collct_init (rcollection)
  
!<description>
  
  ! This routine initialises a collection. It must be applied 
  ! before doing anything to it, just to initialise.
  
!</description>
  
!<inputoutput>
  
  ! The parameter list to initialise.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>
  
!</subroutine>

  ! Allocate memory for the sections. We have at least one unnamed section
  ! and allocate memory for COLLCT_NSECTIONS more named sections in advance.
  ALLOCATE(rcollection%p_Rsections(1+COLLCT_NSECTIONS))
  
  ! Add one section - the unnamed one
  rcollection%isectionCount = 1
  CALL collct_initsection (rcollection%p_Rsections(1), '')
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE collct_done (rcollection)
  
!<description>
  
  ! This routine releases a parameter list. All memory allocated by the
  ! parameter list is released.
  
!</description>
  
!<inputoutput>
  
  ! The parameter list to release.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i

  ! Structure initialised?
  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Warning: trying to release an uninitialised collection.'
    RETURN
  END IF

  ! Loop through the sections and release the content
  DO i=rcollection%isectionCount,1,-1
    CALL collct_donesection(rcollection%p_Rsections(i))
  END DO

  ! Release all section structures
  DEALLOCATE(rcollection%p_Rsections)
  
  ! Mark the structure as 'empty', finish
  rcollection%isectionCount = 0
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION collct_getmaxlevel_direct (rcollection,ssectionName) RESULT (ilevel)

!<description>
  ! Returns the maximum available level in section ssectionname of a collection.
!</description>

!<input>

  ! The collection
  TYPE(t_collection), INTENT(IN) :: rcollection
  
  ! The name of the section
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName
  
!</input>
  
!<result>
  ! The maximum level in rcollection.
!</result>
  
!</function>

  TYPE(t_collctSection), POINTER :: p_rsection
  
  ! Fetch the section number
  CALL collct_fetchsection(rcollection, ssectionName, p_rsection) 
  
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Error: section not found: ',ssectionName
    CALL sys_halt()
  END IF
  
  ilevel = p_rsection%ilevelCount

  END FUNCTION

  ! ***************************************************************************

!<function>

  INTEGER FUNCTION collct_getmaxlevel_indir (rsection) RESULT (ilevel)

!<description>
  ! Returns the maximum available level in section rsection.
  ! A return value of 0 indicates that only level-independent information
  ! is available in the block at the moment.
!</description>

!<input>

  ! The collection
  TYPE(t_collctsection), INTENT(IN) :: rsection
  
!</input>
  
!<result>
  ! The maximum level in rcollection.
!</result>
  
!</function>

  ilevel = rsection%ilevelCount

  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_addlevel_indir (rsection, ilevelid)
  
!<description>
  ! Adds a new level to a section rsection.
  ! ilevelid returns the id of the new level and is >= 1.
!</description>

!<inputoutput>
  
  ! The section where to add the level
  TYPE(t_collctSection), INTENT(INOUT) :: rsection
  
!</inputoutput>
  
!<output>
  
  ! OPTIONAL: The id of the newly added level.
  INTEGER, INTENT(OUT), OPTIONAL :: ilevelid
  
!</output>

!</subroutine>

  ! Add a new level - reallocate the level list if necessary
  IF (.NOT. ASSOCIATED(rsection%p_Rlevels)) THEN
    ALLOCATE(rsection%p_Rlevels(COLLCT_NLEVELS))
  ELSE IF (rsection%ilevelCount .EQ. SIZE(rsection%p_Rlevels)) THEN
    CALL collct_reallocsection (rsection, SIZE(rsection%p_Rlevels)+COLLCT_NLEVELS)
  END IF
  rsection%ilevelCount = rsection%ilevelCount + 1
  
  ! Initialise the new level.
  CALL collct_initlevel(rsection%p_Rlevels(rsection%ilevelCount))

  IF (PRESENT(ilevelid)) ilevelid = rsection%ilevelCount

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_addlevel_direct (rcollection, ilevelid, ssectionName)
  
!<description>
  ! Adds a new level to a section with the name ssectionName.
  ! If ssectionName='', the new level is added to the unnamed section.
  ! ilevelid returns the id of the new level and is >= 1.
!</description>

!<inputoutput>
  
  ! The collection where to add the level
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>
  
!<input>

  ! OPTIONAL: The section name where to add the level. If ='' or not
  ! given, the level is added to the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>  
  
!<output>
  
  ! The id of the newly added level.
  INTEGER, INTENT(OUT) :: ilevelid
  
!</output>

!</subroutine>

  ! local variables
  TYPE(t_collctSection), POINTER :: p_rsection
  
  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    CALL sys_halt()
  END IF
  
  ! Get the section
  IF (PRESENT(ssectionName)) THEN
    CALL collct_fetchsection(rcollection, ssectionName, p_rsection) 
    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      PRINT *,'Error: section not found: ',ssectionName
      CALL sys_halt()
    END IF
  ELSE
    ! unnamed section
    p_rsection => rcollection%p_Rsections(1)
  END IF

  ! Add the level
  CALL collct_addlevel_indir (p_rsection, ilevelid)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_addlevel_all (rcollection)
  
!<description>
  ! Adds a new level to all sections in the collection.
!</description>

!<inputoutput>
  
  ! The collection where to add the level
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i,ilevelid
  
  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    CALL sys_halt()
  END IF
  
  ! Go through all sections
  DO i=1,rcollection%isectionCount
    ! Add the level there
    CALL collct_addlevel_indir (rcollection%p_Rsections(i), ilevelid)  
  END DO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_addsection (rcollection, ssectionName)
  
!<description>
  ! Adds a new section to a collection. The section gets the name ssectionName.
!</description>

!<inputoutput>
  
  ! The parameter list where to add the section.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
  
  ! The name of the section. Must not be ''!
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName
  
!</input>
  
!</subroutine>

  IF (ssectionName .EQ. '') THEN
    PRINT *,'collct_addSection: Empty section name!'
    CALL sys_halt()
  END IF

  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    CALL sys_halt()
  END IF

  ! Add a new section - reallocate the memory if necessary
  IF (.NOT. ASSOCIATED(rcollection%p_Rsections)) THEN
    ALLOCATE(rcollection%p_Rsections(COLLCT_NSECTIONS))
  ELSE IF (rcollection%isectionCount .EQ. SIZE(rcollection%p_Rsections)) THEN
    CALL collct_realloccollection (rcollection, SIZE(rcollection%p_Rsections)+COLLCT_NSECTIONS)
  END IF
  rcollection%isectionCount = rcollection%isectionCount + 1
  
  ! Initialise the new section.
  CALL collct_initsection(rcollection%p_Rsections(rcollection%isectionCount),ssectionName)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_deletesection (rcollection, ssectionName)
  
!<description>
  ! Removes the section ssectionName with all its content from the collection.
  !
  ! Warning: Removing a section from the collection makes all pointers to
  ! sections invalid!
!</description>

!<inputoutput>
  
  ! The parameter list where to add the section.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
  
  ! The name of the section. Must not be ''!
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName
  
!</input>
  
!</subroutine>

  TYPE(t_collctSection), POINTER :: p_rsection
  INTEGER :: isectionIndex,i

  ! Some basic checks
  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    CALL sys_halt()
  END IF

  CALL collct_fetchsection(rcollection, ssectionName, p_rsection, isectionIndex) 
  
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'collct_clearsection: Section does not exist: ',ssectionName
    RETURN
  END IF
  
  ! Remove the section data
  CALL collct_donesection (p_rsection)
  
  ! Now 'relocate' all sections. This makes all pointers to sections
  ! in the main application invalid, but it's rather unlikely that the
  ! application maintains a section pointer in a stage where sections are
  ! modified - hopefully :-)
  
  DO i=isectionIndex+1,rcollection%isectionCount
    rcollection%p_Rsections(i-1) = rcollection%p_Rsections(i)
  END DO
  
  ! Reduce the number of sections in the collection, finish.
  rcollection%isectionCount = rcollection%isectionCount - 1
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_deletevalue (rcollection, sparameter, ilevel, ssectionName) 

!<description>
  ! Deletes the stored parameter 'sparameter' from the collection.
  ! This simply removes a pointer or value from the collection, no memory
  ! is released.
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName
!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  TYPE(t_collctLevel), POINTER :: p_rlevel
  INTEGER :: ilv
  
  ilv = 0
  IF (PRESENT(ilevel)) ilv = ilevel
  
  ! Get the parameter
  IF (PRESENT(ssectionName)) THEN
    CALL collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                       sparameter, p_rvalue) 
  ELSE
    CALL collct_fetchparameter_direct (rcollection, '', ilv, &
                                       sparameter, p_rvalue) 
  END IF
  
  ! Does the value exist?
  IF (ASSOCIATED(p_rvalue)) THEN
  
    ! Ok, we have the parameter.
    ! Deleting it is simple: Clear the name and assign all values
    ! to default values. Then the node is 'dead' and can be reused
    ! later.

    p_rvalue%sname = ''
    p_rvalue%csvalue = " "
    NULLIFY(p_rvalue%p_svalue)
    p_rvalue%ivalue = 0
    p_rvalue%dvalue = 0.0_DP
    NULLIFY(p_rvalue%p_rdiscretisation)
    NULLIFY(p_rvalue%p_rtriangulation)
    NULLIFY(p_rvalue%p_rvectorScalar)
    NULLIFY(p_rvalue%p_rmatrixScalar)
    NULLIFY(p_rvalue%p_rvector)
    NULLIFY(p_rvalue%p_rmatrix)
    NULLIFY(p_rvalue%p_rparlist)
    NULLIFY(p_rvalue%p_rlinearSolver)
    NULLIFY(p_rvalue%p_rdiscreteBC)
    NULLIFY(p_rvalue%p_rboundary)
    NULLIFY(p_rvalue%p_rboundaryConditions)
    NULLIFY(p_rvalue%p_rparser)
    NULLIFY(p_rvalue%p_RfilterChain)
    NULLIFY(p_rvalue%p_rcollection)
    NULLIFY(p_rvalue%p_rhadapt)
    NULLIFY(p_rvalue%p_rafcstab)
    
    IF (ASSOCIATED(p_rvalue%p_Iarray)) DEALLOCATE(p_rvalue%p_Iarray)
    IF (ASSOCIATED(p_rvalue%p_Darray)) DEALLOCATE(p_rvalue%p_Darray)

    p_rvalue%itype = COLLCT_UNDEFINED
    
    ! Modify the level info:
    IF (PRESENT(ssectionName)) THEN
      CALL collct_fetchlevel(rcollection, ssectionName, ilv, p_rlevel)     
    ELSE
      CALL collct_fetchlevel(rcollection, '', ilv, p_rlevel)      
    END IF
    
    ! Decrement the value counter
    p_rLevel%ivalueCount = p_rLevel%ivalueCount - 1
    
    ! Indicate that there's at least one free position now
    p_rLevel%bisFull = .FALSE.
    
  END IF
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_printStatistics (rcollection)

!<description>
  ! This routine prints the current content of the whole collection in
  ! user-readable form to screen.
!</description>  
  
!<input>
  ! The parameter list.
  TYPE(t_collection), INTENT(IN) :: rcollection
!</input>

!</subroutine>

  INTEGER :: isection, ilevel
  
  ! Loop through all sections in the collection
  DO isection=1,rcollection%isectionCount
  
    ! Print the segment name:
    CALL output_line ('['//TRIM(rcollection%p_Rsections(isection)%ssectionName)//']')
    
    ! Print the content of level 0 of the current section
    CALL printlevel(0,rcollection%p_Rsections(isection)%rlevel0)
    
    ! Loop through all levels and print them
    DO ilevel = 1,rcollection%p_Rsections(isection)%ilevelCount
      CALL printlevel(ilevel,rcollection%p_Rsections(isection)%p_Rlevels(ilevel))
    END DO
    
  END DO
  
  CONTAINS
  
    !****************************************************************
    ! Print out the content of a level to screen
    SUBROUTINE printlevel (ilevel,rlevel)
    
    ! Number of the level; 0=level 0
    INTEGER, INTENT(IN) :: ilevel
    
    ! The level structure
    TYPE(t_collctLevel), INTENT(IN) :: rlevel
    
    ! local variables
    INTEGER :: ivalue
    TYPE(t_collctValue), POINTER :: p_rvalue
    
      ! Is there data in the level?
      IF (rlevel%ivalueCount .NE. 0) THEN
      
        ! A level > 0 has a name:
        IF (ilevel .GT. 0) THEN
          CALL output_line ('  [Level '//TRIM(sys_siL(ilevel,8))//']')
        END IF
        
        ! Loop through all values
        DO ivalue = 1,SIZE(rlevel%p_Rvalues)
          ! Get the value structure
          p_rvalue => rlevel%p_Rvalues(ivalue)
          
          ! Is there data in it?
          IF (p_rvalue%itype .NE. COLLCT_UNDEFINED) THEN
          
            ! Print that there is data:
            CALL output_line ('  '//TRIM(p_rvalue%sname)//' (Type: '//&
                              TRIM(sys_siL(p_rvalue%itype,8))//')')
          
          END IF
        END DO
      
      END IF
      
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION collct_queryvalue_indir (rsection, sparameter, ilevel) &
               RESULT (exists)
          
!<description>
  ! Checks whether a parameter sparameter exists in in section rsection.
!</description>  
  
!<result>
  ! >0, if the parameter exists in the section rsection (it's the index of the
  !     parameter within the section),
  ! =0, if the parameter does not exist whithin the section.
!</result>

!<input>
    
  ! The parameter list.
  TYPE(t_collctSection), INTENT(IN), TARGET :: rsection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel
  
!</input>
  
!</function>

  ! local variables
  INTEGER :: ilv
  TYPE(t_collctLevel), POINTER :: p_rlevel
  
  exists = 0
  
  IF (.NOT. PRESENT(ilevel)) THEN
    ilv = 0
  ELSE
    ilv = MAX(ilevel,0)
  END IF
  
  ! Get the right level structure
  IF (ilv .EQ. 0) THEN
    p_rlevel => rsection%rlevel0
  ELSE IF (ilevel .GT. rsection%ilevelCount) THEN
    PRINT *,'collct_queryvalue_indir: Level out of bounds!'
    CALL sys_halt()
  ELSE
    p_rlevel => rsection%p_Rlevels(ilv)
  END IF
  
  ! Search for the parameter in the level
  CALL collct_fetchparameter_indir(p_rlevel, sparameter, exists) 

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION collct_queryvalue_direct (rcollection, sparameter, &
                        ilevel, ssectionName) RESULT (exists)
          
!<description>
  ! Checks whether a parameter sparameter exists in in section.
  ! ilevel and ssectionName are optional parameters.
  ! If ilevel=0 or not given, the parameter is searched in the level-independent
  ! parameter set, otherwise on level ilevel.
  ! If ssectionName='' or not given, the parameter is searched in the unnamed
  ! section, otherwise in section ssectionName.
!</description>  
  
!<result>
  ! >0, if the parameter exists in the section (it's the index of the
  !     parameter within the section),
  ! =0, if the parameter does not exist whithin the section.
!</result>

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(IN) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</function>

  ! local variables
  INTEGER :: ilv
  TYPE(t_collctSection), POINTER :: p_rsection
  
  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    CALL sys_halt()
  END IF

  ! Fall back to the 'indirect' search.
  IF (.NOT. PRESENT(ilevel)) THEN
    ilv = 0
  ELSE
    ilv = MAX(0,ilevel)
  END IF
  
  IF (.NOT. PRESENT(ssectionName)) THEN
    p_rsection => rcollection%p_Rsections(1)
  ELSE
    CALL collct_fetchsection(rcollection, ssectionName, p_rsection) 
    IF (.NOT. ASSOCIATED(p_rsection)) THEN 
      PRINT *,'collct_queryvalue_direct: Section not found'
      CALL sys_halt()
    END IF
  END IF
    
  exists = collct_queryvalue_indir (p_rsection, sparameter, ilv)

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION collct_gettype (rcollection, sparameter, &
                        ilevel, ssectionName, bexists) RESULT (itype)
!<description>
  ! Returns the data type of the parameter sparameter.
!</description>  
  
!<result>
  ! One of the COLLCT_XXXX constants (COLLCT_UNDEFINED,...).
!</result>

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(IN) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  INTEGER :: ilv
  
  ilv = 0
  IF (PRESENT(ilevel)) ilv = ilevel
  
  ! Get the parameter
  IF (PRESENT(ssectionName)) THEN
    CALL collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                       sparameter, p_rvalue) 
  ELSE
    CALL collct_fetchparameter_direct (rcollection, '', ilv, &
                                       sparameter, p_rvalue) 
  END IF
  
  ! Return whether or not that thing exists
  IF (PRESENT(bexists)) bexists = ASSOCIATED(p_rvalue)
  
  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    itype = p_rvalue%itype
  ELSE
    itype = COLLCT_UNDEFINED
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  SUBROUTINE collct_getvalue_struc (rcollection, sparameter, itype, &
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
  ! If bexists is present, it's set to TRUE/FALSE to indicate whether
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! Expected type for that parameter.
  ! =COLLCT_UNDEFINED: suppress type check.
  INTEGER, INTENT(IN) :: itype
  
  ! Whether to add the parameter if it does not exists.
  ! =TRUE: add the parameter if necessary.
  ! =FALSE: throw an error that the parameter does not exist.
  LOGICAL, INTENT(IN) :: badd
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists
  
  ! A pointer to the structure with name sparameter.
  ! Points to NULL() if the value does not exist.
  TYPE(t_collctValue), POINTER :: p_rvalue

!</output>

!</function>

  ! local variables
  INTEGER :: ilv
  
  ilv = 0
  IF (PRESENT(ilevel)) ilv = ilevel
  
  ! Get the parameter
  IF (PRESENT(ssectionName)) THEN
    CALL collct_fetchparameter_direct (rcollection, ssectionName, ilv, &
                                       sparameter, p_rvalue) 
  ELSE
    CALL collct_fetchparameter_direct (rcollection, '', ilv, &
                                       sparameter, p_rvalue) 
  END IF
  
  ! Return whether or not that thing exists.
  ! If it does not exist, create it or throw an error.
  IF (PRESENT(bexists)) THEN
    bexists = ASSOCIATED(p_rvalue)
  ELSE
    IF ((.NOT. badd) .AND. (.NOT. ASSOCIATED(p_rvalue))) THEN
      PRINT *,'collct_getvalue: Parameter '//TRIM(sparameter)//' not found!'
      CALL sys_halt()
    END IF
  END IF

  IF (ASSOCIATED(p_rvalue)) THEN

    ! Type check
    IF (itype .NE. COLLCT_UNDEFINED) THEN
      ! Throw an error if the type is wrong. Otherwise, get the value.
      IF (p_rvalue%itype .NE. itype) THEN
        PRINT *,'collct_getvalue: Wrong type!'
        PRINT *,'Parameter: '//TRIM(sparameter)// ' is of type ',p_rvalue%itype,&
                ' but expected as ',itype,'!'
        CALL sys_halt()
      END IF
    END IF

  ELSE
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              itype, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              itype, ilv, p_rvalue)
      END IF
      IF (PRESENT(bexists)) bexists = .TRUE.
    END IF    
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<function>

  CHARACTER FUNCTION collct_getvalue_char (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) RESULT (value)
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_CHARACTER,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value = p_rvalue%csvalue
  ELSE
    value = " "
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_getvalue_string (rcollection, sparameter, value, &
                                     ilevel, ssectionName, bexists)
!<description>
  ! Returns the the parameter sparameter as string.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! The string.
  CHARACTER(LEN=*) :: value

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists
  
!</output>

!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  INTEGER :: i,j

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_STRING,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    i = MIN(SIZE(p_rvalue%p_svalue),LEN(value))
    value = ''
    DO j=1,i
      value(j:j) = p_rvalue%p_svalue(j)
    END DO
  ELSE
    value = ''
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<function>

  INTEGER(I32) FUNCTION collct_getvalue_int (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) RESULT(value)
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGER,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value = p_rvalue%ivalue
  ELSE
    value = 0
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_getvalue_intarr (rcollection, sparameter, value, &
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
  INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: value
!</inputoutput>

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>
  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists
!</output>

!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGERARR,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    IF (SIZE(p_rvalue%p_Iarray) .NE. SIZE(value)) THEN
      CALL output_line ('Destination array has the wrong length!',&
          OU_CLASS_ERROR,OU_MODE_STD,'collct_getvalue_intarr')
      CALL sys_halt()
    END IF
    value = p_rvalue%p_Iarray
  ELSE
    value = 0
  END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<function>

  REAL(DP) FUNCTION collct_getvalue_real (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) RESULT(value)
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_REAL,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value = p_rvalue%dvalue
  ELSE
    value = 0
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  SUBROUTINE collct_getvalue_realarr (rcollection, sparameter, value, &
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
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: value
!</inputoutput>

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_REALARR,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    IF (SIZE(p_rvalue%p_Darray) .NE. SIZE(value)) THEN
      CALL output_line ('Destination array has the wrong length!',&
          OU_CLASS_ERROR,OU_MODE_STD,'collct_getvalue_realarr')
      CALL sys_halt()
    END IF
    value = p_rvalue%p_Darray
  ELSE
    value = 0
  END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_discr (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a discretisation 
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_spatialDiscretisation), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCR,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rdiscretisation
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_bldiscr (rcollection, sparameter, &
                                    ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a block discretisation 
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_blockDiscretisation), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BLDISCR,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rblockDiscretisation
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_tria (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a triangulation 
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_triangulation), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_TRIA,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rtriangulation
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_domain (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a triangulation 
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_boundary), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARY,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rboundary
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_bc (rcollection, sparameter, &
                               ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a triangulation 
  ! structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_boundaryConditions), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARYCOND,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rboundaryConditions
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_vecsca (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar vector.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_vectorScalar), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAVECTOR,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rvectorScalar
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_matsca (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar matrix.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_matrixScalar), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAMATRIX,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rmatrixScalar
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_vec (rcollection, sparameter, &
                                ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a block vector.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_vectorBlock), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKVECTOR,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rvector
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_mat (rcollection, sparameter, &
                               ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar matrix.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_matrixBlock), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKMATRIX,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rmatrix
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_linsol (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a linear solver
  ! configuration node.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_linsolNode), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_LINSOL,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rlinearSolver
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_discbc (rcollection, sparameter, &
                                  ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a discrete boundary
  ! condition structure
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_discreteBC), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCRBC,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rdiscreteBC
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_plst (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a parameter list.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_parlist), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_PARAMETERS,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rparlist
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_ilvpsc (rcollection, sparameter, &
                                   ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a scalar interlevel
  ! projection structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_interlevelProjectionScalar), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJSC,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rilvProjectionSc
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_ilvp (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a (block) interlevel
  ! projection structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_interlevelProjectionBlock), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJ,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rilvProjection
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_coll (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a collection object.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_collection), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_COLLECTION,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rcollection
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_pars (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a parser object.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_fparser), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_FPARSER,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rparser
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_geom (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a geometry object.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_geometryObject), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_GEOMETRY,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rgeometry
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_fchn (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a filter chain.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_filterChain), DIMENSION(:), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_FILTERCHAIN,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_RfilterChain
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_hadapt (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a grid adaptation structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_hadapt), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_HADAPT,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rhadapt
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  FUNCTION collct_getvalue_afcstab (rcollection, sparameter, &
                                 ilevel, ssectionName, bexists) RESULT(value)
!<description>
  ! Returns the the parameter sparameter as pointer to a stabilisation structure.
  ! An error is thrown if the value is of the wrong type.
!</description>  
  
!<result>
  ! The value of the parameter.
  ! A standard value if the value does not exist.
!</result>
  
  TYPE(t_afcstab), POINTER :: value

!<input>
    
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!<output>

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists

!</output>

!</function>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue

  ! Get the pointer to the parameter
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_HADAPT,&
                              .FALSE.,p_rvalue, ilevel, bexists, ssectionName)

  ! Return the quantity
  IF (ASSOCIATED(p_rvalue)) THEN
    value => p_rvalue%p_rafcstab
  ELSE
    NULLIFY(value)
  END IF

  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_char (rcollection, sparameter, value, badd, &
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  CHARACTER, INTENT(IN) :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_CHARACTER,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! Set the value
  p_rvalue%csvalue = value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_string (rcollection, sparameter, value, badd, &
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  CHARACTER(LEN=*), INTENT(IN) :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists
  INTEGER :: i

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_STRING,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! Set the value; slightly more complicated for strings
  ALLOCATE(p_rvalue%p_svalue(LEN(value)))
  DO i=1,LEN(value)
    p_rvalue%p_svalue(i) = value (i:i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_int (rcollection, sparameter, value, badd, &
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  INTEGER, INTENT(IN) :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGER,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! Set the value
  p_rvalue%ivalue = value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_intarr (rcollection, sparameter, value, badd, &
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTEGERARR,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! (Re-)allocate memory for the value if necessary.
  IF (.NOT. ASSOCIATED(p_rvalue%p_Iarray)) THEN
    ALLOCATE(p_rvalue%p_Iarray(SIZE(value)))
  ELSE IF (SIZE(p_rvalue%p_Iarray) .NE. SIZE(value)) THEN
    DEALLOCATE(p_rvalue%p_Iarray)
    ALLOCATE(p_rvalue%p_Iarray(SIZE(value)))
  END IF
  
  ! Set the value
  p_rvalue%p_Iarray = value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_real (rcollection, sparameter, value, badd, &
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  REAL(DP), INTENT(IN) :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_REAL,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! Set the value
  p_rvalue%dvalue = value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_realarr (rcollection, sparameter, value, badd, &
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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  REAL(DP), DIMENSION(:),INTENT(IN) :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_REALARR,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! (Re-)allocate memory for the value if necessary.
  IF (.NOT. ASSOCIATED(p_rvalue%p_Darray)) THEN
    ALLOCATE(p_rvalue%p_Darray(SIZE(value)))
  ELSE IF (SIZE(p_rvalue%p_Darray) .NE. SIZE(value)) THEN
    DEALLOCATE(p_rvalue%p_Darray)
    ALLOCATE(p_rvalue%p_Darray(SIZE(value)))
  END IF
  
  ! Set the value
  p_rvalue%p_Darray = value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_discr (rcollection, sparameter, value, badd, &
                                    ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCR,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! Set the value
  p_rvalue%p_rdiscretisation => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_bldiscr (rcollection, sparameter, value, badd, &
                                      ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BLDISCR,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! Set the value
  p_rvalue%p_rblockDiscretisation => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_tria (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_triangulation), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_TRIA,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)

  ! Set the value
  p_rvalue%p_rtriangulation => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_domain (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_boundary), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARY,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rboundary => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_bc (rcollection, sparameter, value, badd, &
                                 ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_boundaryConditions), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BOUNDARYCOND,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rboundaryConditions => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_vecsca (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_vectorScalar), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAVECTOR,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rvectorScalar => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_matsca (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_matrixScalar), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_SCAMATRIX,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rmatrixScalar => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_vec (rcollection, sparameter, value, badd, &
                                  ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKVECTOR,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rvector => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_mat (rcollection, sparameter, value, badd, &
                                  ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_BLKMATRIX,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rmatrix => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_parlst (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_parlist), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_PARAMETERS,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rparlist => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_linsol (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_linsolNode), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_LINSOL,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rlinearSolver => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_discbc (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_discreteBC), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_DISCRBC,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rdiscreteBC => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_ilvpsc (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_interlevelProjectionScalar), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJSC,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rilvProjectionSc => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_ilvp (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_interlevelProjectionBlock), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_INTERLVPRJ,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rilvProjection => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_coll (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_collection), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_COLLECTION,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rcollection => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_pars (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_fparser), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_FPARSER,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rparser => value

  END SUBROUTINE

! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_geom (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_geometryObject), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_GEOMETRY,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rgeometry => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_fchn (rcollection, sparameter, value, badd, &
                                   ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_filterChain), DIMENSION(:), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_FILTERCHAIN,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_RfilterChain => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_hadapt (rcollection, sparameter, value, badd, &
                                     ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_hadapt), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_FILTERCHAIN,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rhadapt => value

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE collct_setvalue_afcstab (rcollection, sparameter, value, badd, &
                                      ilevel, ssectionName) 
!<description>
  ! Stores a pointer to 'value' using the parametre name 'sparameter'.
  ! If the parameter does not exist, the behaviour depends on the 
  ! parameter badd:
  !  badd=false: an error is thrown,
  !  badd=true : the parameter is created at the position defined by
  !              ilevel and ssectionName (if given). When the position
  !              defined by these variables does not exist, an error is thrown
!</description>  
  
!<inputoutput>
  
  ! The parameter list.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
  
!</inputoutput>

!<input>
    
  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
  ! The value of the parameter.
  TYPE(t_afcstab), INTENT(IN), TARGET :: value
  
  ! Whether to add the variable if it does not exist.
  ! =false: don't add the variable, throw an error
  ! =true : add the variable
  LOGICAL, INTENT(IN) :: badd

  ! OPTIONAL: The level where to search.
  ! If =0 or not given, the search is in the level-independent parameter block.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel

  ! OPTIONAL: The section name where to search.
  ! If ='' or not given, the search is in the unnamed section.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssectionName

!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  LOGICAL :: bexists

  ! Get the pointer to the parameter. Add the parameter if necessary
  CALL collct_getvalue_struc (rcollection, sparameter, COLLCT_FILTERCHAIN,&
                              badd,p_rvalue, ilevel, bexists, ssectionName)
  
  ! Set the value
  p_rvalue%p_rafcstab => value

  END SUBROUTINE

END MODULE
