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
!# - character string
!# - double real
!# - integer 
!# - a vector - scalar and block type
!# - a matrix - scalar and block type
!# - triangulation structures
!# - discrete BC structures
!# - geometry information structures
!# - linear solver configurations
!# - a parameter list
!# - another collection
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
!# 10.) collct_detelevalue
!#      -> Deletes a value from the list.
!#
!# 11.) collct_printStatistics 
!#      -> Prints out the current content of the structure to the terminal
!#
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
!#      rcollection%Iquickaccess(1) = 1       ! Number of boundary component e.g.
!#      rcollection%Iquickaccess(2) = 3       ! Number of Dirichlet-segment e.g.
!#      rcollection%Dquickaccess(1) = 1.0_DP  ! Inflow speed e.g.
!#      CALL bcasm_discretiseBC (...,callback,rcollection)
!#      ...
!#    END PROGRAM
!#
!#    SUBROUTINE callback(...,rbcRegion,...,rcollection,...,Dvalues)
!#      ...
!#      iboundaryComponent = rcollection%Iquickaccess(1)
!#      iinflowSegment = rcollection%Iquickaccess(2)
!#      dinfowSpeed    = rcollection%Dquickaccess(1)
!#      IF ((iboundaryComponent .EQ. rbcRegion%iboundCompIdx) .AND.
!#          (iinflowSegment .EQ. rbcRegion%iboundSegIdx)) THEN
!#        Dvalues(1) = dinflowSpeed
!#      END IF
!#      ...
!#    END SUBROUTINE
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
!</constantblock>

!<constantblock description='Type identifier for values'>
  
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

  ! Discretisation structure
  INTEGER, PARAMETER :: COLLCT_DISCR     = 5

  ! Triangulation structure
  INTEGER, PARAMETER :: COLLCT_TRIA      = 6

  ! A scalar vector
  INTEGER, PARAMETER :: COLLCT_SCAVECTOR = 7

  ! A scalar matrix
  INTEGER, PARAMETER :: COLLCT_SCAMATRIX = 8

  ! A block vector
  INTEGER, PARAMETER :: COLLCT_BLKVECTOR = 9

  ! A block matrix
  INTEGER, PARAMETER :: COLLCT_BLKMATRIX = 10
  
  ! A parameter structure
  INTEGER, PARAMETER :: COLLCT_PARAMETERS = 11

  ! A linear solver structure
  INTEGER, PARAMETER :: COLLCT_LINSOL     = 12

  ! Discretised-boundary-conditions structure
  INTEGER, PARAMETER :: COLLCT_DISCRBC    = 13

  ! A domain
  INTEGER, PARAMETER :: COLLCT_BOUNDARY   = 14

  ! Scalar analytic boundary conditions
  INTEGER, PARAMETER :: COLLCT_BOUNDARYCOND = 15

  ! The collection structure itself
  INTEGER, PARAMETER :: COLLCT_COLLECTION = 16

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
    TYPE(t_boundary), POINTER      :: p_rdomain => NULL()

    ! Pointer to scalar boundary conditions
    TYPE(t_boundaryConditions), POINTER      :: p_rboundaryConditions => NULL()

    ! Pointer to a collection structure
    TYPE(t_collection), POINTER              :: p_rcollection => NULL()

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
  
    ! Quick-access array for integers. This is a short, temporary
    ! array that can be used by the application to save intermediate
    ! values (e.g. some information that is to pass and to be accessed
    ! by callback routines). 
    ! No long life information should be stored here. Long-life
    ! information should be named and added as such to the collection 
    ! directly. The collection itself does not maintain this array!
    INTEGER(I32), DIMENSION(COLLCT_QALENGTH) :: IquickAccess
  
    ! Quick-access array for reals. This is a short, temporary
    ! array that can be used by the application to save intermediate
    ! values (e.g. some information that is to pass and to be accessed
    ! by callback routines). 
    ! No long life information should be stored here. Long-life
    ! information should be named and added as such to the collection 
    ! directly. The collection itself does not maintain this array!
    REAL(DP), DIMENSION(COLLCT_QALENGTH) :: DquickAccess
  
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
    ALLOCATE(p_Rvalues(SIZE(p_Rvalues)+COLLCT_NVALUES))
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
        STOP
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
  p_rvalue%sname = sparamName
  p_rvalue%itype = itype
  
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
    STOP
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
  ! and return a pointer to the section - or NULL)( if the section does 
  ! not exist.

  SUBROUTINE collct_fetchsection(rcollection, sname, p_rsection) 

  ! The section.
  TYPE(t_collection), INTENT(IN) :: rcollection
  
  ! The parameter name to look for. 
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! A pointer to the section or NULL() if it does not exist.
  TYPE(t_collctSection), POINTER :: p_rsection
  
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
    STOP
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
    STOP
  END IF
  
  ! Get the section
  IF (PRESENT(ssectionName)) THEN
    CALL collct_fetchsection(rcollection, ssectionName, p_rsection) 
    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      PRINT *,'Error: section not found: ',ssectionName
      STOP
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
    STOP
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
    STOP
  END IF

  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    STOP
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

    p_rvalue%itype = COLLCT_UNDEFINED
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
    NULLIFY(p_rvalue%p_rdomain)
    NULLIFY(p_rvalue%p_rboundaryConditions)
    NULLIFY(p_rvalue%p_rcollection)
    
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
    PRINT *,'['//TRIM(rcollection%p_Rsections(isection)%ssectionName)//']'
    
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
          PRINT *,'  [Level ',ilevel,']'
        END IF
        
        ! Loop through all values
        DO ivalue = 1,SIZE(rlevel%p_Rvalues)
          ! Get the value structure
          p_rvalue => rlevel%p_Rvalues(ivalue)
          
          ! Is there data in it?
          IF (p_rvalue%itype .NE. COLLCT_UNDEFINED) THEN
          
            ! Print that there is data:
            PRINT *,'  '//TRIM(p_rvalue%sname)//' (Type: ',&
                    p_rvalue%itype,')'
          
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
  
!</description>
  
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
    STOP
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
  
!</description>
  
!</function>

  ! local variables
  INTEGER :: ilv
  TYPE(t_collctSection), POINTER :: p_rsection
  
  IF (rcollection%isectionCount .EQ. 0) THEN
    PRINT *,'Error: Collection not initalised!'
    STOP
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
      STOP
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_CHARACTER) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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

  ! The string.
  CHARACTER(LEN=*) :: value

  ! OPTIONAL: Returns TRUE if the variable exists, FALSE otherwise.
  ! There's no error thrown if a variable does not exist.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists
  
!</output>

!</subroutine>

  ! local variables
  TYPE(t_collctValue), POINTER :: p_rvalue
  INTEGER :: ilv,i,j
  
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_STRING) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_INTEGER) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
    value = p_rvalue%ivalue
  ELSE
    value = 0
  END IF

  END FUNCTION

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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_REAL) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
    value = p_rvalue%dvalue
  ELSE
    value = 0
  END IF

  END FUNCTION

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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_DISCR) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
    value => p_rvalue%p_rdiscretisation
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_TRIA) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_BOUNDARY) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
    value => p_rvalue%p_rdomain
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_BOUNDARYCOND) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_SCAVECTOR) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_SCAMATRIX) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_BLKVECTOR) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_BLKMATRIX) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_PARAMETERS) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_DISCRBC) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_DISCRBC) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
    value => p_rvalue%p_rparlist
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
    ! Throw an error if the type is wrong. Otherwise, get the value.
    IF (p_rvalue%itype .NE. COLLCT_COLLECTION) THEN
      PRINT *,'Wrong type! Parameter: ',sparameter, ' is of type ',p_rvalue%itype
      STOP
    END IF
    
    value => p_rvalue%p_rcollection
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_CHARACTER, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_CHARACTER, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  INTEGER :: ilv,i
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_STRING, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_STRING, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  INTEGER(I32), INTENT(IN) :: value
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_INTEGER, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_INTEGER, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
  ! Set the value
  p_rvalue%ivalue = value

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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_REAL, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_REAL, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
  ! Set the value
  p_rvalue%dvalue = value

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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_DISCR, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_DISCR, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
  ! Set the value
  p_rvalue%p_rdiscretisation => value

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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_TRIA, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_TRIA, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_BOUNDARY, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_BOUNDARY, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
  ! Set the value
  p_rvalue%p_rdomain => value

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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_BOUNDARYCOND, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_BOUNDARYCOND, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_SCAVECTOR, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_SCAVECTOR, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_SCAMATRIX, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_SCAMATRIX, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_BLKVECTOR, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_BLKVECTOR, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_BLKMATRIX, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_BLKMATRIX, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_PARAMETERS, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_PARAMETERS, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_LINSOL, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_LINSOL, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_DISCRBC, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_DISCRBC, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
  ! Set the value
  p_rvalue%p_rdiscreteBC => value

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
  
  ! Add the value if necessary
  IF (.NOT. ASSOCIATED(p_rvalue)) THEN
    IF (badd) THEN
      IF (PRESENT(ssectionName)) THEN
        CALL collct_addvalue (rcollection, ssectionName, sparameter, &
                              COLLCT_COLLECTION, ilv, p_rvalue)
      ELSE
        CALL collct_addvalue (rcollection, '', sparameter, &
                              COLLCT_COLLECTION, ilv, p_rvalue)
      END IF
    ELSE
      PRINT *,'Error: Parameter ',sparameter,' does not exist!'
      STOP
    END IF
  END IF
  
  ! Set the value
  p_rvalue%p_rcollection => value

  END SUBROUTINE

END MODULE
