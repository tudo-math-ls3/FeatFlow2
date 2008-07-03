!##############################################################################
!# ****************************************************************************
!# <name> paramlist </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a dynamic parameter list in the memory. It allows
!# to read .INI parameter from disc, parses it and saves the content into a
!# structure. The parameters can later be accessed by the program.
!#
!# A typical .INI file has the following structure:
!# -------------------------snip------------------------------
!# # this is a data file
!#
!# parameter1 = data1
!# parameter2 = data2
!#
!# # now a section follows
!#
!# [Section1]
!# parameter3 = 'This is a string of ''multiple'' words'
!# parameter4 = data3
!# 
!# [section2]
!# parameter5 = data4    # This is a comment
!#
!# [section3]
!# parameterlist(4)=     # An array consisting of 4 strings
!#   data-line1
!#   data-line2
!#   data-line3
!#   data-line4
!# -------------------------snip------------------------------
!# the .INI file is build up by different sections.
!# Each section starts with the section name enclosed by
!# brackets ('[...]'), followed by a list of parameters 
!# if the form "name=parameter". There's one unnamed parameter
!# block at the beginning of the file.
!#
!# There's at most one parameter per line. The parameter names 
!# are case-insensitive. Allowed characters for parameter names
!# are: 'A'..'Z', '0'..'9', '-', '_', '(', ')'.
!#
!# The parameter data is always a string consisting of
!# - one word without spaces or
!# - multiple words, enclosed by apostrophes.
!# Spaces in the parameter data are ignored except when they
!# are enclosed by apostrophes. 
!# A parameter name followed by "(n)" is assumed to be an array.
!# The following n nonempty lines define string values for all
!# n entries in the array.
!# 
!# Empty lines are ignored. When there is a '#' character in
!# a line not enclosed by apostrophes, the rest of the line is
!# ignored as a comment.
!#
!# The following routines can be used to maintain a parameter
!# list:
!#  1.) parlst_init 
!#       -> Initialises an empty parameter list
!#
!#  2.) parlst_readfromfile
!#      -> Reads the content of a .INI file into a parameter list.
!#
!#  3.) parlst_clear
!#      -> Cleans up a parameter list, removes all parameters.
!#
!#  4.) parlst_done
!#      -> Cleans up a parameter list, releases all allocated memory
!#         from the heap
!#
!#  5.) parlst_querysection
!#      -> Determines whether or not a section exists
!#
!#  6.) parlst_addsection
!#      -> Adds a new section
!#
!#  7.) parlst_queryvalue
!#      -> Determines whether or not a parameter exists
!#
!#  8.) parlst_querysubstrings
!#      -> Returns the number of substrings in a parameter
!#
!#  9.) parlst_getvalue_string
!#      parlst_getvalue_int
!#      parlst_getvalue_double
!#      -> Get the string/int/real value of a parameter from the parameter list
!#
!# 10.) parlst_addvalue
!#     -> Adds a new parameter to the parameter list
!#
!# 11.) parlst_setvalue
!#      -> Modifies the value of a parameter in the list
!#
!# 12.) parlst_getStringRepresentation
!#      -> Creates a string representation of the parameter list.
!#
!# 13.) parlst_info 
!#      -> Print the parameter list to the terminal
!#
!# </purpose>
!##############################################################################

MODULE paramlist

  USE fsystem
  USE io
  USE genoutput
  
  IMPLICIT NONE

!<constants>

  !<constantblock>
  
  ! Maximum length of a section name.
  INTEGER, PARAMETER :: PARLST_MLSECTION = 64

  ! Maximum length of parameter names: 32 characters
  INTEGER, PARAMETER :: PARLST_MLNAME = 32

  ! Maximum length of parameter data: 256 characters
  INTEGER, PARAMETER :: PARLST_MLDATA = 256
  
  ! Minimum number of free parameter 'slots' per parameter section.
  ! If there are too many parameters in a parameter section, the
  ! structure is dynamically extended in terms of PARLST_NPARSPERBLOCK
  ! entries.
  INTEGER, PARAMETER :: PARLST_NPARSPERBLOCK = 32

  ! Minimum number of parameter sections.
  ! If there are too many parameter sections in a parameter block, the
  ! structure is dynamically extended in terms of PARLST_NSECTIONS
  ! entries.
  INTEGER, PARAMETER :: PARLST_NSECTIONS = 8

  ! Maximum length of a line in a INI file. Lines longer than this
  ! are truncated.
  INTEGER, PARAMETER :: PARLST_LENLINEBUF = 1024
  
  ! Comment character
  CHARACTER, PARAMETER :: PARLST_COMMENT = "#"

  !</constantblock>

!</constants>

!<types>
  
  !<typeblock>
  
  ! This structure realises a value associated to a parameter name.
  ! A value consists of one string or an array of strings.
  
  TYPE t_parlstValue
  
    PRIVATE
    
    ! Number of strings. If set to 0, the value consists of one
    ! string, to be found in svalue. If > 0, there are nsize
    ! strings to be found in p_Sentry.
    INTEGER :: nsize = 0
    
    ! Single string; contains the value in case nsize=0
    CHARACTER(LEN=PARLST_MLDATA) :: sentry = ''
    
    ! Array of strings in case nsize>0
    CHARACTER(LEN=PARLST_MLDATA), DIMENSION(:), POINTER :: p_Sentry => NULL()
  
  END TYPE
  
  !</typeblock>
  
  !<typeblock>
  
  ! This structure realises a parameter section. It contains an
  ! array with parameter names and an array with parameter values
  ! to these names. The arrays are dynamically allocated. 
  
  TYPE t_parlstSection
  
    PRIVATE
  
    ! The name of the section.
    CHARACTER(LEN=PARLST_MLSECTION) :: ssectionName = ''
    
    ! Actual number of parameters in this section.
    INTEGER :: iparamCount = 0
    
    ! A list of parameter names. Each name contains PARLST_MLNAME
    ! characters.
    CHARACTER(LEN=PARLST_MLNAME), DIMENSION(:), POINTER :: p_Sparameters => NULL()
    
    ! A list of t_parlstValue structures corresponding to the parameters
    ! in p_Sparameters.
    TYPE(t_parlstValue), DIMENSION(:), POINTER :: p_Rvalues
    
  END TYPE
  
  !</typeblock>
  
  !<typeblock>
  
  ! This structure realises a parameter list. Parameters can be read into
  ! it from a file. Parameters can be obtained from the structure using
  ! the query/get routines.
  
  TYPE t_parlist
  
    PRIVATE
  
    ! Actual number of sections in the parameter list. There's at least
    ! one section - the unnamed section. If this value is =0, the parameter
    ! list is not initialised.
    INTEGER :: isectionCount = 0
    
    ! A list of sections. The first section is always the unnamed section.
    TYPE(t_parlstSection), DIMENSION(:), POINTER :: p_Rsections => NULL()
    
  END TYPE
  
  !</typeblock>

!</types>

  PRIVATE :: parlst_initsection, parlst_reallocsection, parlst_realloclist
  PRIVATE :: parlst_fetchparameter,parlst_readlinefromfile,parlst_parseline
  
  INTERFACE parlst_queryvalue
    MODULE PROCEDURE parlst_queryvalue_direct
    MODULE PROCEDURE parlst_queryvalue_indir
  END INTERFACE

  INTERFACE parlst_querysubstrings
    MODULE PROCEDURE parlst_querysubstrings_direct
    MODULE PROCEDURE parlst_querysubstrings_indir
  END INTERFACE

  INTERFACE parlst_addvalue
    MODULE PROCEDURE parlst_addvalue_direct
    MODULE PROCEDURE parlst_addvalue_indir
  END INTERFACE

  INTERFACE parlst_setvalue
    MODULE PROCEDURE parlst_setvalue_fetch
    MODULE PROCEDURE parlst_setvalue_indir
    MODULE PROCEDURE parlst_setvalue_direct
  END INTERFACE

  INTERFACE parlst_getvalue_string
    MODULE PROCEDURE parlst_getvalue_string_fetch
    MODULE PROCEDURE parlst_getvalue_string_indir
    MODULE PROCEDURE parlst_getvalue_string_direct
  END INTERFACE

  INTERFACE parlst_getvalue_int
    MODULE PROCEDURE parlst_getvalue_int_fetch
    MODULE PROCEDURE parlst_getvalue_int_indir
    MODULE PROCEDURE parlst_getvalue_int_direct
  END INTERFACE

  INTERFACE parlst_getvalue_double
    MODULE PROCEDURE parlst_getvalue_double_fetch
    MODULE PROCEDURE parlst_getvalue_double_indir
    MODULE PROCEDURE parlst_getvalue_double_direct
  END INTERFACE

CONTAINS
  
  ! ***************************************************************************

  ! Internal subroutine: Initialise a newly created parameter section.
  
  SUBROUTINE parlst_initsection (rparlstSection,sname)
  
  TYPE(t_parlstSection), INTENT(INOUT) :: rparlstSection
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! Simply allocate the pointers with an empty list
  ALLOCATE(rparlstSection%p_Sparameters(PARLST_NPARSPERBLOCK))
  ALLOCATE(rparlstSection%p_Rvalues(PARLST_NPARSPERBLOCK))
  
  ! and set the section name
  rparlstSection%ssectionName = sname
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Reallocate a section.
  ! This increases the size of a parameter section by reallocation of the
  ! arrays.
  
  SUBROUTINE parlst_reallocsection (rparlstSection, inewsize)
  
  ! The section to reallocate.
  TYPE(t_parlstSection), INTENT(INOUT) :: rparlstSection
  
  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the section should be able to handle.
  INTEGER, INTENT(IN) :: inewsize
  
  ! local variables
  
  INTEGER :: sz,oldsize

  ! Pointers to new lists for replacing the old.
  CHARACTER(LEN=PARLST_MLNAME), DIMENSION(:), POINTER :: p_Sparameters 
  TYPE(t_parlstValue), DIMENSION(:), POINTER :: p_Rvalues
  
  oldsize = SIZE(rparlstSection%p_Sparameters)
  sz = MAX(oldsize,inewsize)
  
  IF (SIZE(rparlstSection%p_Sparameters) .EQ. sz) RETURN ! nothing to do
  
  ! Allocate the pointers for the new lists
  ALLOCATE(p_Sparameters(sz))
  ALLOCATE(p_Rvalues(sz))
  
  ! Copy the content of the old ones
  p_Sparameters(1:oldsize) = rparlstSection%p_Sparameters (1:oldsize)
  p_Rvalues(1:oldsize) = rparlstSection%p_Rvalues (1:oldsize)
  
  ! Throw away the old arrays, replace by the new ones
  DEALLOCATE(rparlstSection%p_Rvalues)
  DEALLOCATE(rparlstSection%p_Sparameters)
  
  rparlstSection%p_Sparameters => p_Sparameters
  rparlstSection%p_Rvalues => p_Rvalues
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Release a section.
  ! Removes all temporary memory that is allocated by a section.
  
  SUBROUTINE parlst_releasesection (rparlstSection)
  
  ! The section to release.
  TYPE(t_parlstSection), INTENT(INOUT) :: rparlstSection
  
  ! local variables
  INTEGER :: i
  
  ! Loop through all values in the current section if there is
  ! an array-value. Release them.
  DO i=SIZE(rparlstSection%p_Rvalues),1,-1
    IF (rparlstSection%p_Rvalues(i)%nsize .GT. 0) THEN
      DEALLOCATE(rparlstSection%p_Rvalues(i)%p_Sentry)
    END IF
  END DO
  
  ! Remove the content of the section.
  DEALLOCATE(rparlstSection%p_Rvalues)
  DEALLOCATE(rparlstSection%p_Sparameters)
  rparlstSection%iparamCount = 0
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Reallocate the section list
  ! This increases the size of a section list by reallocation of the
  ! arrays.
  
  SUBROUTINE parlst_realloclist (rparlist, inewsize)
  
  ! The section list to reallocate.
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
  
  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the section should be able to handle.
  INTEGER, INTENT(IN) :: inewsize
  
  ! local variables
  
  INTEGER :: sz

  ! Pointers to new lists for replacing the old.
  TYPE(t_parlstSection), DIMENSION(:), POINTER :: p_Rsections
  
  ! Allocate the pointers for the new lists
  ALLOCATE(p_Rsections(inewsize))
  
  sz = MIN(SIZE(rparlist%p_Rsections),inewsize)
  
  ! Copy the content of the old ones
  p_Rsections(1:sz) = rparlist%p_Rsections (1:sz)
  
  ! Throw away the old arrays, replace by the new ones
  DEALLOCATE(rparlist%p_Rsections)
  
  rparlist%p_Rsections => p_Rsections
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
  ! Internal subroutine: Convert a character to upper case.
  
  SUBROUTINE parlst_toupper (str) 
  
  ! The string that is to make uppercase
  CHARACTER(LEN=*), INTENT(INOUT) :: str
  
!  CHARACTER(LEN=26), PARAMETER :: lowc = 'abcdefghijklmnopqrstuvwxyz'
!  CHARACTER(LEN=26), PARAMETER :: upc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!  
!  INTEGER :: i
!  
!  i = INDEX(lowc,c)
!  IF (i .NE. 0) THEN
!    c = upc(i:i)
!  END IF

  INTEGER, PARAMETER :: up2low = IACHAR("a") - IACHAR("A")
  INTEGER :: i
  CHARACTER    :: c
  
  DO i=1,LEN(str)
    c = str(i:i)
    IF ((c .GE. "a") .AND. (c .LE. "z")) THEN
      str(i:i) = ACHAR (IACHAR(c) - up2low)
    END IF
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
  ! Internal subroutine: Search in a section for a parameter
  ! and return the index - or 0 if the parameter does not exist.

  SUBROUTINE parlst_fetchparameter(rsection, sname, iparamnum) 

  ! The section.
  TYPE(t_parlstSection), INTENT(IN) :: rsection
  
  ! The parameter name to look for. Must be uppercase.
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! The number of the parameter in the list or 0 if it does not exist.
  INTEGER, INTENT(OUT) :: iparamnum
  
  ! local variables
  INTEGER :: i
  
  iparamnum = 0
  
  ! If the parameter list is empty, the section does not exist for sure
  IF (rsection%iparamCount .EQ. 0) RETURN
  
  ! Loop through all sections to see if the section exists
  DO i=1,rsection%iparamCount
    IF (rsection%p_Sparameters(i) .EQ. sname) THEN
      iparamnum = i
      RETURN
    END IF
  END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE parlst_init (rparlist)
  
!<description>
  
  ! This routine initialises a parameter list. It must be applied to a
  ! parameter list structure before doing anything to it, just to initialise.
  
!</description>
  
!<inputoutput>
  
  ! The parameter list to initialise.
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
  
!</inputoutput>
  
!</subroutine>

  ! Set the section-count to 1.
  rparlist%isectionCount = 1
  
  ! Allocate a first set of sections
  ALLOCATE(rparlist%p_Rsections(PARLST_NSECTIONS))
  
  ! Initialise the first section - it's the unnamed one.
  CALL parlst_initsection (rparlist%p_Rsections(1),'')
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE parlst_clear (rparlist)
  
!<description>
  ! This routine cleans up a parameter list. All parameters in rparlist are
  ! removed.
!</description>
  
!<inputoutput>
  ! The parameter list to clean up.
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
!</inputoutput>
  
!</subroutine>

    ! Clean up = done+reinit. We make that simple here...
    CALL parlst_done (rparlist)
    CALL parlst_init (rparlist)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE parlst_done (rparlist)
  
!<description>
  
  ! This routine releases a parameter list. All memory allocated by the
  ! parameter list is released.
  
!</description>
  
!<inputoutput>
  
  ! The parameter list to release.
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i

  ! Probably nothing to do
  IF (rparlist%isectionCount .EQ. 0) RETURN

  ! Loop through the parameter lists and release the content
  DO i=rparlist%isectionCount,1,-1
    CALL parlst_releasesection (rparlist%p_Rsections(i))
  END DO

  ! Release all sections
  DEALLOCATE(rparlist%p_Rsections)
  
  ! Mark the structure as 'empty', finish
  rparlist%isectionCount = 0

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE parlst_querysection(rparlist, sname, p_rsection) 

!<description>

  ! Searches for a section and return a pointer to it -
  ! or NULL() of the section does not exist.
  
!</description>

!<input>

  ! The parameter list to scan for the section.
  TYPE(t_parlist), INTENT(IN) :: rparlist
  
  ! The section name to look for. 
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
!</input>
  
!<output>
  
  ! A pointer to the section.
  TYPE(t_parlstSection), POINTER :: p_rsection
  
!</output>
  
!</subroutine>
  
  ! local variables
  INTEGER :: i
  CHARACTER(LEN=PARLST_MLSECTION) :: sectionname
  
  NULLIFY(p_rsection)
  
  ! If the parameter list is empty, the section does not exist for sure
  IF (rparlist%isectionCount .EQ. 0) RETURN
  
  ! If the section name is '', return a pointer to the first section.
  IF (sname .EQ. '') THEN
    p_rsection => rparlist%p_Rsections(1)
    RETURN
  END IF
  
  ! Create the upper-case section name
  sectionname = ADJUSTL(sname)
  CALL parlst_toupper (sectionname)

  ! Loop through all sections to see if the section exists
  DO i=1,rparlist%isectionCount
    IF (rparlist%p_Rsections(i)%ssectionName .EQ. sectionname) THEN
      p_rsection => rparlist%p_Rsections(i)
      RETURN
    END IF
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE parlst_addsection (rparlist, sname)
  
!<description>
  
  ! Adds a section with the name sname to the list of sections in the
  ! parameter list rparlist. The name must NOT contain brackets ('[',']')
  ! in front and at the end!
  
!</description>

!<inputoutput>
  
  ! The parameter list where to add the section.
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
  
!</inputoutput>

!<input>
  
  ! The section name to add - without brackets in front and at the end!
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
!</input>
  
!</subroutine>

  ! local variables
  CHARACTER(LEN=PARLST_MLSECTION) :: sectionname
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    CALL sys_halt()
  END IF
  
  ! Create the upper-case section name
  sectionname = ADJUSTL(sname)
  CALL parlst_toupper (sectionname)
  
  ! Add a new section - reallocate the section list if necessary
  IF (rparlist%isectionCount .EQ. SIZE(rparlist%p_Rsections)) THEN
    CALL parlst_realloclist (rparlist, SIZE(rparlist%p_Rsections)+PARLST_NSECTIONS)
  END IF
  rparlist%isectionCount = rparlist%isectionCount + 1
  
  ! Initialise the new section.
  CALL parlst_initsection(rparlist%p_Rsections(rparlist%isectionCount),sectionname)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION parlst_queryvalue_indir (rsection, sparameter) &
               RESULT (exists)
          
!<description>
  ! Checks whether a parameter sparameter exists in the section rsection.
!</description>
  
!<result>
  ! The index of the parameter in the section ssection or =0, if the
  ! parameter does not exist within the section.
!</result>

!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  exists = 0
  
  IF (sparameter .EQ. '') THEN
    PRINT *,'Empty parameter name!'
    CALL sys_halt()
  END IF
  
  ! Create the upper-case parameter name
  paramname = ADJUSTL(sparameter)
  CALL parlst_toupper (paramname)
  
  ! Get the parameter index into 'exists', finish.
  CALL parlst_fetchparameter(rsection, paramname, exists)
  
  END FUNCTION
  
  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION parlst_queryvalue_direct (rparlist, ssectionName, sparameter) &
               RESULT (exists)
          
!<description>
  ! Checks whether a parameter sparameter exists in the section ssectionname
  ! in the parameter list rparlist.
!</description>
  
!<result>
  ! The index of the parameter in the section ssectionName or =0, if the
  ! parameter does not exist within the section.
!</result>

!<input>
    
  ! The parameter list.
  TYPE(t_parlist), INTENT(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  exists = 0
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    CALL sys_halt()
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Section not found'
    RETURN
  END IF
  
  ! Search for the parameter
  exists = parlst_queryvalue_indir (p_rsection, sparameter)

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION parlst_querysubstrings_indir (rsection, sparameter) &
               RESULT (iresult)
          
!<description>
  ! Returns the number of substrings of a parameter.
!</description>
  
!<result>
  ! The number of substrings of parameter sparameter in section rsection.
!</result>

!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  INTEGER :: idx
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  IF (sparameter .EQ. '') THEN
    PRINT *,'Empty parameter name!'
    CALL sys_halt()
  END IF
  
  ! Create the upper-case parameter name
  paramname = ADJUSTL(sparameter)
  CALL parlst_toupper (paramname)
  
  ! Get the parameter index into 'idx', finish.
  CALL parlst_fetchparameter(rsection, paramname, idx)
  
  ! Return number of substrings
  iresult = rsection%p_Rvalues(idx)%nsize
  
  END FUNCTION
  
  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION parlst_querysubstrings_direct (rparlist, ssectionName, sparameter) &
               RESULT (iresult)
          
!<description>
  ! Checks whether a parameter sparameter exists in the section ssectionname
  ! in the parameter list rparlist.
!</description>
  
!<result>
  ! The index of the parameter in the section ssectionName or =0, if the
  ! parameter does not exist within the section.
!</result>

!<input>
    
  ! The parameter list.
  TYPE(t_parlist), INTENT(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The parameter name to search for.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  INTEGER :: idx
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    CALL sys_halt()
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Section not found'
    RETURN
  END IF
  
  ! Get the parameter index
  idx = parlst_queryvalue_indir (p_rsection, sparameter)

  ! Return number of substrings
  iresult = p_rsection%p_Rvalues(idx)%nsize

  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_string_indir (rsection, &
                                           sparameter, svalue, sdefault, &
                                           isubstring)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, sdefault is returned.
  ! If sdefault is not given, an error will be thrown.
  !
  ! If the value is an array of strings, the optional parameter isubstring>=0
  ! allows to specify the number of the substring to be returned; 
  ! isubstring=0 returns the value directly
  ! behind the '=' sign in the line of the parameter, isubstring>0 returns
  ! the array-entry in the lines below the parameter.
  !
  ! When ommitting isubstring, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! Optional: A default value
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sdefault
  
  ! Optional: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  INTEGER, INTENT(IN), OPTIONAL :: isubstring
  
!</input>
  
!<output>

  ! The value of the parametzer
  CHARACTER(LEN=*), INTENT(OUT) :: svalue
  
!</output>

!</subroutine>

  ! local variables
  INTEGER :: i,isub
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  IF (sparameter .EQ. '') THEN
    PRINT *,'Empty parameter name!'
    CALL sys_halt()
  END IF

  ! Create the upper-case parameter name
  paramname = ADJUSTL(sparameter)
  CALL parlst_toupper (paramname)
  
  ! Get the parameter index into 'exists', finish.
  CALL parlst_fetchparameter(rsection, paramname, i)
  
  IF (i .EQ. 0) THEN
    IF (PRESENT(sdefault)) THEN
      svalue = sdefault
    ELSE
      PRINT *,'Parameter ',TRIM(paramname),' does not exist!'
      CALL sys_halt()
    END IF
  ELSE
    ! Depending on isubstring, return either the 'headline' or one
    ! of the substrings.
    isub = 0
    IF (PRESENT(isubstring)) isub = isubstring
  
    IF ((isub .LE. 0) .OR. (isub .GT. rsection%p_Rvalues(i)%nsize)) THEN
      svalue = rsection%p_Rvalues(i)%sentry
    ELSE
      svalue = rsection%p_Rvalues(i)%p_Sentry(isub)
    END IF
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_string_fetch (rsection, &
                                           iparameter, svalue, bexists,&
                                           isubstring)
!<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  ! If bexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If bexists is given, it will be set to TRUE if the parameter number
  ! iparameter exists, otherwise it will be set to FALSE and svalue=''.
  !  
  ! If the value is an array of strings, the optional parameter isubstring>=0
  ! allows to specify the number of the substring to be returned; 
  ! isubstring=0 returns the value directly
  ! behind the '=' sign in the line of the parameter, isubstring>0 returns
  ! the array-entry in the lines below the parameter.
  !
  ! When ommitting isubstring, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The number of the parameter.
  INTEGER, INTENT(IN) :: iparameter

  ! Optional: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  INTEGER, INTENT(IN), OPTIONAL :: isubstring

!</input>
  
!<output>

  ! The value of the parameter
  CHARACTER(LEN=*), INTENT(OUT) :: svalue
  
  ! Optional: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists
  
!</output>

!</subroutine>

  INTEGER :: isub

  ! Check if iparameter is out of bounds. If yes, probably
  ! throw an error.
  
  IF ((iparameter .LT. 0) .OR. (iparameter .GT. rsection%iparamCount)) THEN
  
    IF (.NOT. PRESENT(bexists)) THEN 
      PRINT *,'Error. Parameter ',iparameter,' does not exist!'
      CALL sys_halt()
    ELSE
      svalue = ''
      bexists = .FALSE.
      RETURN
    END IF
  
  END IF
  
  ! Get the parameter value.
  ! Depending on isubstring, return either the 'headline' or one
  ! of the substrings.
  isub = 0
  IF (PRESENT(isubstring)) isub = isubstring

  IF ((isub .LE. 0) .OR. &
      (isub .GT. rsection%p_Rvalues(iparameter)%nsize)) THEN
    svalue = rsection%p_Rvalues(iparameter)%sentry
  ELSE
    svalue = rsection%p_Rvalues(iparameter)%p_Sentry(isub)
  END IF
  
  IF (PRESENT(bexists)) bexists = .TRUE.

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_string_direct (rparlist, ssectionName, &
                                            sparameter, svalue, sdefault,&
                                            isubstring)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, sdefault is returned.
  ! If sdefault is not given, an error will be thrown.
  !
  ! If the value is an array of strings, the optional parameter isubstring>=0
  ! allows to specify the number of the substring to be returned; 
  ! isubstring=0 returns the value directly
  ! behind the '=' sign in the line of the parameter, isubstring>0 returns
  ! the array-entry in the lines below the parameter.
  !
  ! When ommitting isubstring, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The parameter list.
  TYPE(t_parlist), INTENT(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! Optional: A default value
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sdefault
  
  ! Optional: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  INTEGER, INTENT(IN), OPTIONAL :: isubstring

!</input>
  
!<output>

  ! The value of the parametzer
  CHARACTER(LEN=*), INTENT(OUT) :: svalue
  
!</output>

!</subroutine>

  ! local variables
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    CALL sys_halt()
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    IF (PRESENT(sdefault)) THEN
      svalue = sdefault
      RETURN
    ELSE
      PRINT *,'Section not found'
      CALL sys_halt()
    END IF
  END IF

  ! Get the value
  CALL parlst_getvalue_string_indir (p_rsection, sparameter, svalue, sdefault,&
                                     isubstring)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_double_indir (rsection, sparameter, dvalue, ddefault)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, idefault is returned.
  ! If idefault is not given, an error will be thrown.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! Optional: A default value
  REAL(DP), INTENT(IN), OPTIONAL :: ddefault
  
!</input>
  
!<output>

  ! The value of the parametzer
  REAL(DP), INTENT(OUT) :: dvalue
  
!</output>

!</subroutine>

  ! local variables
  CHARACTER (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  IF (PRESENT(ddefault)) THEN
    WRITE (sdefault,*) ddefault
    CALL parlst_getvalue_string_indir (rsection, sparameter, svalue, sdefault)
  ELSE
    CALL parlst_getvalue_string_indir (rsection, sparameter, svalue)
  END IF
  
  READ(svalue,*) dvalue

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE parlst_getvalue_double_fetch (rsection, iparameter, dvalue, bexists)

!<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  ! If bexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If bexists is given, it will be set to TRUE if the parameter number
  ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The number of the parameter.
  INTEGER, INTENT(IN) :: iparameter

!</input>
  
!<output>

  ! The value of the parameter
  INTEGER, INTENT(OUT) :: dvalue
  
  ! Optional: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists
  
!</output>

!</subroutine>

  ! local variables
  CHARACTER (LEN=PARLST_MLDATA) :: svalue
  
  svalue = '0.0E0'
  CALL parlst_getvalue_string_fetch (rsection, &
                                     iparameter, svalue, bexists)
  READ(svalue,*) dvalue

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_double_direct (rparlist, ssectionName, &
                                            sparameter, dvalue, ddefault)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, ddefault is returned.
  ! If ddefault is not given, an error will be thrown.
  
!</description>
  
!<input>
    
  ! The parameter list.
  TYPE(t_parlist), INTENT(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! Optional: A default value
  REAL(DP), INTENT(IN), OPTIONAL :: ddefault
  
!</input>
  
!<output>

  ! The value of the parametzer
  REAL(DP), INTENT(OUT) :: dvalue
  
!</output>

!</subroutine>

  ! local variables
  CHARACTER (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  IF (PRESENT(ddefault)) THEN
    WRITE (sdefault,*) ddefault
    CALL parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, sdefault)
  ELSE
    CALL parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue)
  END IF
  
  READ(svalue,*) dvalue

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_int_indir (rsection, sparameter, ivalue, idefault)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, idefault is returned.
  ! If idefault is not given, an error will be thrown.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! Optional: A default value
 INTEGER, INTENT(IN), OPTIONAL :: idefault
  
!</input>
  
!<output>

  ! The value of the parametzer
  INTEGER, INTENT(OUT) :: ivalue
  
!</output>

!</subroutine>

  ! local variables
  CHARACTER (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  IF (PRESENT(idefault)) THEN
    WRITE (sdefault,*) idefault
    CALL parlst_getvalue_string_indir (rsection, sparameter, svalue, sdefault)
  ELSE
    CALL parlst_getvalue_string_indir (rsection, sparameter, svalue)
  END IF
  
  READ(svalue,*) ivalue

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_int_fetch (rsection, iparameter, ivalue, bexists)
!<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  ! If bexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If bexists is given, it will be set to TRUE if the parameter number
  ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The number of the parameter.
  INTEGER, INTENT(IN) :: iparameter

!</input>
  
!<output>

  ! The value of the parameter
  INTEGER, INTENT(OUT) :: ivalue
  
  ! Optional: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  LOGICAL, INTENT(OUT), OPTIONAL :: bexists
  
!</output>

!</subroutine>

  ! local variables
  CHARACTER (LEN=PARLST_MLDATA) :: svalue
  
  svalue = '0'
  CALL parlst_getvalue_string_fetch (rsection, &
                                     iparameter, svalue, bexists)
  READ(svalue,*) ivalue

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_int_direct (rparlist, ssectionName, &
                                         sparameter, ivalue, idefault)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, idefault is returned.
  ! If idefault is not given, an error will be thrown.
  
!</description>
  
!<input>
    
  ! The parameter list.
  TYPE(t_parlist), INTENT(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! Optional: A default value
  INTEGER, INTENT(IN), OPTIONAL :: idefault
  
!</input>
  
!<output>

  ! The value of the parametzer
  INTEGER, INTENT(OUT) :: ivalue
  
!</output>

!</subroutine>

  ! local variables
  CHARACTER (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  IF (PRESENT(idefault)) THEN
    WRITE (sdefault,*) idefault
    CALL parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, sdefault)
  ELSE
    CALL parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue)
  END IF
  
  READ(svalue,*) ivalue

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE parlst_addvalue_indir (rsection, sparameter, svalue, nsubstrings)
  
!<description>
  ! Adds a parameter to a section rsection.
!</description>
  
!<inputoutput> 
    
  ! The section where to arr the parameter
  TYPE(t_parlstSection), INTENT(INOUT) :: rsection
  
!</inputoutput>

!<input>

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! The value of the parameter
  CHARACTER(LEN=*), INTENT(IN) :: svalue
  
  ! Optional: Number of substrings. This allows a parameter to have
  ! multiple substrings, which can be accessed via the 'isubstring'
  ! parameter in the GET-routines.
  INTEGER, INTENT(IN), OPTIONAL :: nsubstrings
  
!</input>

!</subroutine>

  ! local variables
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  ! Create the upper-case parameter name
  paramname = ADJUSTL(sparameter)
  CALL parlst_toupper (paramname)

  ! Enough space free? Otherwise reallocate the parameter list
  IF (rsection%iparamCount .EQ. SIZE(rsection%p_Sparameters)) THEN
    CALL parlst_reallocsection (rsection, SIZE(rsection%p_Sparameters)+PARLST_NPARSPERBLOCK)
  END IF

  ! Add the parameter - without any adjustment of the 'value' string
  rsection%iparamCount = rsection%iparamCount + 1  
  
  rsection%p_Sparameters(rsection%iparamCount) = paramname
  rsection%p_Rvalues(rsection%iparamCount)%sentry = svalue
  
  ! Add a list for the substrings if the parameter should have substrings.
  IF (PRESENT(nsubstrings)) THEN
    IF (nsubstrings .GT. 0) THEN
      ALLOCATE(rsection%p_Rvalues(rsection%iparamCount)%p_Sentry(nsubstrings))
      rsection%p_Rvalues(rsection%iparamCount)%nsize = nsubstrings
    END IF
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_addvalue_direct (rparlist, ssectionName, sparameter, svalue,&
                                     nsubstrings)
!<description>
  
  ! Adds a parameter to a section with name ssectionName in the parameter list
  ! rparlist. If ssectionName='', the parameter is added to the unnamed
  ! section.
  
!</description>
  
!<inputoutput> 
    
  ! The parameter list.
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
  
!</inputoutput>

!<input>

  ! The section name - '' identifies the unnamed section.
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! The value of the parameter
  CHARACTER(LEN=*), INTENT(IN) :: svalue
  
  ! Optional: Number of substrings. This allows a parameter to have
  ! multiple substrings, which can be accessed via the 'isubstring'
  ! parameter in the GET-routines.
  INTEGER, INTENT(IN), OPTIONAL :: nsubstrings

!</input>

!</subroutine>

  ! local variables
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    CALL sys_halt()
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Section not found'
    RETURN
  END IF

  ! Add the parameter 
  
  CALL parlst_addvalue_indir (p_rsection, sparameter, svalue, nsubstrings)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE parlst_setvalue_fetch (rsection, iparameter, svalue, iexists,&
                                    isubstring)
  
!<description>
  
  ! Modifies the value of a parameter in the section rsection.
  ! The value of parameter iparameter in the section rsection is modified.
  ! If iexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If iexists is given, it will be set to YES if the parameter number
  ! iparameter exists and was modified, otherwise it will be set to NO.
  !
  ! isubstring allows to specify the numer of a substring of the parameter to
  ! change. If ommitted or = 0, the 'headline' directly behind the '='
  ! sign of the line 'name=value' is modified. Otherwise, the corresponding
  ! substring is changed.
  
!</description>
  
!<inputoutput> 
    
  ! The section where to arr the parameter
  TYPE(t_parlstSection), INTENT(INOUT) :: rsection
  
!</inputoutput>

!<input>

  ! The parameter name.
  INTEGER, INTENT(IN) :: iparameter

  ! The new value of the parameter
  CHARACTER(LEN=*), INTENT(IN) :: svalue
  
  ! Optional: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  INTEGER, INTENT(IN), OPTIONAL :: isubstring

!</input>

!<output>

  ! Optional parameter. Is set to YES/NO, depending on whether
  ! the parameter exists.
  INTEGER, INTENT(OUT), OPTIONAL :: iexists

!</output>

!</subroutine>

  INTEGER :: isub

  ! Check if iparameter is out of bounds. If yes, probably
  ! throw an error.
  
  IF ((iparameter .LT. 0) .OR. (iparameter .GT. rsection%iparamCount)) THEN
  
    IF (.NOT. PRESENT(iexists)) THEN 
      PRINT *,'Error. Parameter ',iparameter,' does not exist!'
      CALL sys_halt()
    ELSE
      iexists = NO
      RETURN
    END IF
  
  END IF

  ! Depending on isubstring, change either the 'headline' or one
  ! of the substrings.
  isub = 0
  IF (PRESENT(isubstring)) isub = isubstring

  IF ((isub .LE. 0) .OR. &
      (isub .GT. rsection%p_Rvalues(iparameter)%nsize)) THEN
    rsection%p_Rvalues(iparameter)%sentry = svalue
  ELSE
    rsection%p_Rvalues(iparameter)%p_Sentry(isub) = svalue
  END IF

  IF (PRESENT(iexists)) iexists = YES

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE parlst_setvalue_indir (rsection, sparameter, svalue, isubstring)
  
!<description>
  
  ! Modifies the value of a parameter in the section rsection.
  ! If the parameter does not exist, an error is thrown.
  !
  ! isubstring allows to specify the numer of a substring of the parameter to
  ! change. If ommitted or = 0, the 'headline' directly behind the '='
  ! sign of the line 'name=value' is modified. Otherwise, the corresponding
  ! substring is changed.
  
!</description>
  
!<inputoutput> 
    
  ! The section where to arr the parameter
  TYPE(t_parlstSection), INTENT(INOUT) :: rsection
  
!</inputoutput>

!<input>

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! The new value of the parameter
  CHARACTER(LEN=*), INTENT(IN) :: svalue
  
  ! Optional: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  INTEGER, INTENT(IN), OPTIONAL :: isubstring

!</input>

!</subroutine>

  ! local variables
  INTEGER :: i,isub
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  ! Create the upper-case parameter name
  paramname = ADJUSTL(sparameter)
  CALL parlst_toupper (paramname)

  ! Get the parameter position
  i = parlst_queryvalue_indir (rsection, paramname)
  
  IF (i .EQ. 0) THEN
    PRINT *,'Parameter ',paramname,' does not exist, cannot be modified!'
    CALL sys_halt()
  ELSE 
  
    ! Depending on isubstring, change either the 'headline' or one
    ! of the substrings.
    isub = 0
    IF (PRESENT(isubstring)) isub = isubstring

    IF ((isub .LE. 0) .OR. (isub .GT. rsection%p_Rvalues(i)%nsize)) THEN
      rsection%p_Rvalues(i)%sentry = svalue
    ELSE
      rsection%p_Rvalues(i)%p_Sentry(isub) = svalue
    END IF
  
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_setvalue_direct (rparlist, ssectionName, sparameter, svalue,&
                                     isubstring)
!<description>
  
  ! Modifies the value of a parameter in the section with name ssectionName
  ! in the parameter list rparlist.
  ! If the parameter does not exist, an error is thrown.
  !
  ! isubstring allows to specify the numer of a substring of the parameter to
  ! change. If ommitted or = 0, the 'headline' directly behind the '='
  ! sign of the line 'name=value' is modified. Otherwise, the corresponding
  ! substring is changed.
  
!</description>
  
!<inputoutput> 
    
  ! The parameter list.
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
  
!</inputoutput>

!<input>

  ! The section name - '' identifies the unnamed section.
  CHARACTER(LEN=*), INTENT(IN) :: ssectionName

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! The new value of the parameter
  CHARACTER(LEN=*), INTENT(IN) :: svalue
  
  ! Optional: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  INTEGER, INTENT(IN), OPTIONAL :: isubstring

!</input>

!</subroutine>

  ! local variables
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    CALL sys_halt()
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Section not found'
    RETURN
  END IF

  ! Set the parameter 
  
  CALL parlst_setvalue_indir (p_rsection, sparameter, svalue, isubstring)

  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Read a line from a text file.
  
  SUBROUTINE parlst_readlinefromfile (iunit, sdata, ilinelen, ios)
  
  ! The unit where to read from; must be connected to a file.
  INTEGER, INTENT(IN) :: iunit
  
  ! The string where to write data to
  CHARACTER(LEN=*), INTENT(OUT) :: sdata
  
  ! Length of the output
  INTEGER, INTENT(OUT) :: ilinelen
  
  ! Status of the reading process. Set to a value <> 0 if the end
  ! of the file is reached.
  INTEGER, INTENT(OUT) :: ios
  
  ! local variables
  INTEGER :: eol
  CHARACTER :: c
  
  sdata = ''
  ilinelen = 0
  
  ! Read the data - as long as the line/file does not end.
  eol = NO
  ios = 0
  DO WHILE ((ios .EQ. 0) .AND. (eol .EQ. NO))
    
    ! Read a character.
    ! Unfortunately, Fortran forces me to use this dirty GOTO
    ! to decide processor-independently whether the line or
    ! the record ends.
    READ (unit=iunit,fmt='(A1)',iostat=ios,advance='NO', end=10, eor=20) c
    GOTO 30
    
10  CONTINUE
    ! End of file. 
    ios = -1
    GOTO 30
    
20  CONTINUE
    ! End of record = END OF LINE.
    eol = YES

    ! Set error flag back to 0.
    ios = 0
    
30  CONTINUE    
    ! Don't do anything in case of an error
    IF (ios .EQ. 0) THEN
    
      ilinelen = ilinelen + 1
      sdata (ilinelen:ilinelen) = c
    
    END IF
  
  END DO
  
  END SUBROUTINE

  ! ***************************************************************************

  ! Internal subroutine: Parse a text line.
  ! This parses the text line sdata.
  ! Return values:
  !  ityp = 0 -> The line is a comment
  !  ityp = 1 -> The line is a section. ssecname is the name of the section
  !              without '[]'
  !  ityp = 2 -> The line is a parameter. sparamname is the uppercase
  !              parameter name. svalue is the value of the parameter,
  !              trimmed and left adjusted.
  !  ityp = 3 -> Line is the beginning of a multi-valued parameter.
  !              The next isubstring lines contain additional substrings.
  !  ityp = 4 -> Line is a substring of a multi-valued parameter.
  
  SUBROUTINE parlst_parseline (sdata, ityp, isubstring, ilinenum, &
                               ssecname, sparamname, svalue)
  
  ! The line to be parsed
  CHARACTER(LEN=*), INTENT(IN) :: sdata
  
  ! The typ of the line
  INTEGER, INTENT(OUT) :: ityp

  ! input: =0: parse line as parameter. isubstring is changed to a value > 0
  !            is the parameter has multiple values attached.
  !        >0: parse line as substring of a multi-valued parameter, not 
  !            containing a leading 'name='.
  ! output: If the 'headline' of a multi-valued parameter is read, isubstring is
  !         changed to the number of substrings (the k in 'name(k)=...').
  !         Otherwise unchanged.
  INTEGER, INTENT(INOUT) :: isubstring

  ! Line number
  INTEGER, INTENT(IN) :: ilinenum
  
  ! Section name, if it's a section
  CHARACTER(LEN=*), INTENT(INOUT) :: ssecname
  
  ! Parameter name, if it's a parameter
  CHARACTER(LEN=*), INTENT(INOUT) :: sparamname
  
  ! Parameter value, if it's a parameter
  CHARACTER(LEN=*), INTENT(INOUT) :: svalue
  
  ! local variables
  INTEGER :: i,j1,j2,ltr
  CHARACTER(LEN=PARLST_LENLINEBUF) :: sbuf,slen
  
    ityp = 0
    
    ! Do we have data in sdata?
    IF (sdata .EQ. '') RETURN
    
    ! Copy the input string - left adjusted - and get the string length
    sbuf = ADJUSTL(sdata)
    
    ! Should we parse the line as first line of a parameter or as substring
    ! of a multi-valued parameter?
    IF (isubstring .EQ. 0) THEN
    
      ! Standard parameter or section header.
      !    
      ! Do we start with '[' and end with ']'?
      IF (sbuf(1:1) .EQ. "[") THEN
      
        ! Find the final ']'.
        DO ltr = 1,LEN(sbuf)
          IF (sbuf(ltr:ltr) .EQ. "]") EXIT
        END DO
        
        IF (sbuf(ltr:ltr) .NE. ']') THEN
          PRINT *,'Wrong syntax of section name. Line ',ilinenum,':'
          PRINT *,sbuf
          CALL sys_halt()
        END IF
        
        ! Get the section name
        ssecname = sbuf(2:ltr-1)
        ityp = 1
        RETURN
        
      ELSE IF (sbuf(1:1) .EQ. PARLST_COMMENT) THEN
      
        ! Comment sign
        RETURN
        
      ELSE
      
        ! Must be a parameter. Get the length of the string without comment
        ! at the end.
        CALL linelength(sbuf, ltr)
        
        ! ltr=0 means: empty line. Ignore that.
        IF (ltr .EQ. 0) RETURN
        
        ! Is there a '(..)' that is indicating a multi-valued parameter?
        j1 = INDEX(sbuf(1:ltr),'(')
        j2 = INDEX(sbuf(1:ltr),')')

        ! Is there a '=' sign?
        i = INDEX(sbuf(1:ltr),'=')

        IF (i .EQ. 0) THEN
          PRINT *,'Invalid parameter syntax. Line ',ilinenum,':'
          CALL sys_halt()
        END IF
      
        IF ((j1 .EQ. 0) .OR. (j2 .LE. j1)) THEN
        
          ityp = 2
          
          ! Get the name of the parameter
          sparamname = ADJUSTL(sbuf(1:i-1))
          
          ! Get the parameter value
          svalue = ADJUSTL(sbuf(i+1:ltr))
          
        ELSE
        
          ! Probably multi-valued parameter with substrings in the 
          ! following lines.

          ! Get the name of the parameter
          sparamname = ADJUSTL(sbuf(1:j1-1))

          ! Get the parameter value
          svalue = ADJUSTL(sbuf(i+1:ltr))

          ! Get the length of the parameter list.
          slen = sbuf (j1+1:MIN(j2-1,LEN(slen)))

          isubstring = 0
          READ(slen,*) isubstring
          
          IF (isubstring .LE. 0) THEN
            ! Oh, only one line. User want's to cheat :-)
            isubstring = 0
            
            ityp = 2
          ELSE
            ! Real multi-valued parameter.
            ityp = 3
          END IF
        
        END IF
      
      END IF
      
    ELSE
      
      ! Substring of a multi-valued parameter.
      IF (sbuf(1:1) .EQ. PARLST_COMMENT) THEN
      
        ! Comment sign
        RETURN
        
      ELSE
       
        ! Must be a parameter. Get the length of the string without comment
        ! at the end.
        CALL linelength(sbuf, ltr)
        
        ! ltr=0 means: empty line. Ignore that.
        IF (ltr .EQ. 0) RETURN
        
        ityp = 4
        
        ! Get the parameter value. Don't get a parameter name; there is none.
        svalue = ADJUSTL(sbuf(1:ltr))
       
      END IF
    
    END IF
  
  CONTAINS
    
    ! Sub-subroutine: find the length of the line, removing comments
    ! at the end.
    
    SUBROUTINE linelength (sdata, l)
    
    ! The string to parse. Must not be ''!
    CHARACTER(LEN=*), INTENT(IN) :: sdata
    
    ! The index of the last character without any comment at the end.
    INTEGER, INTENT(OUT) :: l
    
    ! local variables
    LOGICAL :: bflag   ! Set to true if we are in apostroph mode
    INTEGER :: lsdata
    
    bflag = .FALSE.
    
    ! Go through all characters
    l = 0
    lsdata = LEN(sdata)
    DO WHILE (l .LT. lsdata)
      
      ! next character
      l = l+1

      ! A comment character while we are not in apostroph mode? Stop.
      IF ((.NOT. bflag) .AND. (sdata(l:l) .EQ. PARLST_COMMENT)) THEN
        l = l-1
        EXIT
      END IF
      
      ! An apostroph? 
      IF (sdata(l:l) .EQ. "'") THEN
      
        ! Switch the apostroph mode.
        ! Btw.: Two subsequent apostrophes will switch the mode off and on again.
        bflag = .NOT. bflag
      
      END IF

    END DO
    
    END SUBROUTINE
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE parlst_readfromfile (rparlist, sfilename)
  
!<description>
  
  ! This routine parses a text file for data of the INI-file form.
  ! sfilename must be the name of a file on the hard disc.
  ! The parameters read from the file are added to the parameter list
  ! rparlist, which has to be initialised with parlst_init before
  ! calling the routine.
  ! Remark: When adding parameters/sections to rparlist, the routine
  !   does not check whether the parameters/sections already exist.
  !   Adding a parameter/section which exists does not result in an error -
  !   the second instance of the parameter/section will simply become
  !   inaccessible, since the GET-routines always return the first
  !   instance!
  
!</description>

!<inputoutput> 
    
  ! The parameter list which is filled with data from the file
  TYPE(t_parlist), INTENT(INOUT) :: rparlist
  
!</inputoutput>
  
!<input>
  
  ! The filename of the file to read.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
  
!</input>
  
!</subroutine>

  ! local variables
  INTEGER :: iunit,ios,isbuflen,ityp,ilinenum,isubstring,nsubstrings,iparpos
  TYPE(t_parlstSection), POINTER :: p_currentsection
  CHARACTER(LEN=PARLST_LENLINEBUF) :: sdata
  CHARACTER(LEN=PARLST_MLSECTION) :: ssectionname
  CHARACTER(LEN=PARLST_MLNAME) :: sparname
  CHARACTER(LEN=PARLST_MLDATA) :: svalue
  
  ! Try to open the file
  CALL io_openFileForReading(sfilename, iunit)
  
  ! Oops...
  IF (iunit .EQ. -1) THEN
    PRINT *,'Error opening .INI file.' 
    CALL sys_halt()
  END IF
  
  ! Start adding parameters to the unnamed section
  p_currentsection => rparlist%p_Rsections(1)
  
  ! Read all lines from the file
  ios = 0
  ilinenum = 0
  isubstring = 0
  nsubstrings = 0
  DO WHILE (ios .EQ. 0) 
    
    ! Read a line from the file into sbuf
    CALL parlst_readlinefromfile (iunit, sdata, isbuflen, ios)
    ilinenum = ilinenum + 1
    
    IF (isbuflen .NE. 0) THEN
    
      ! Parse the line
      CALL parlst_parseline (sdata, ityp, nsubstrings, ilinenum, ssectionname, &
                             sparname, svalue)  
      
      SELECT CASE (ityp)
      CASE (1)
        ! A new section name. Add a section, set the current section
        ! to the new one.
        CALL parlst_addsection (rparlist, ssectionname)
        p_currentsection => rparlist%p_Rsections(rparlist%isectionCount)
        
      CASE (2)
        ! A new parameter. Add it to the current section.
        CALL parlst_addvalue (p_currentsection, sparname, svalue)
        
      CASE (3)
        ! 'Headline' of a multi-valued parameter. Add the parameter with
        ! isubstring subvalues
        CALL parlst_addvalue (p_currentsection, sparname, svalue, nsubstrings)
        
        ! Fetch the parameter for later adding of subvalues.
        iparpos = parlst_queryvalue(p_currentsection, sparname)
        
        ! isubstring counts the current readed substring.
        ! Set it to 0, it will be increased up to nsubstrings in 'case 4'.
        isubstring = 0
        
      CASE (4)
        ! Increase number of current substring
        isubstring = isubstring + 1
        
        ! Sub-parameter of a multi-valued parameter. Add the value to
        ! the last parameter that was added in case 3.
        CALL parlst_setvalue_fetch (p_currentsection, iparpos, svalue, &
                                    isubstring=isubstring)
                                    
        ! Decrement the substring counter. If we reach 0, parlst_parseline
        ! continues to parse standard parameters.
        nsubstrings = nsubstrings - 1
        
      ! Other cases: comment.
      END SELECT
    
    END IF
  
  END DO
  
  ! Close the file, finish.
  CLOSE (iunit)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE parlst_getStringRepresentation (rparlist, p_sconfiguration)
  
!<description>
  ! Creates a string representation of the given parameter list rparlist.
  ! p_sconfiguration will be created as "array[1..*] of char" on
  ! the heap containing this representation. The memory must be manually 
  ! released by the caller using DEALLOCATE when finished using the string
  ! representation.
!</description>
  
!<input> 
  ! The parameter list which is filled with data from the file
  TYPE(t_parlist), INTENT(IN) :: rparlist
!</input>

!<output>
  ! A pointer to a character array containing all lines of the parameter list.
  ! Points to NULL() if there is no data in the parameter list.
  ! Each line is terminated by NEWLINE.
  ! If there is data, a new pointer is allocated for this on the heap. 
  ! The user must manually release the memory when finished using it.
  CHARACTER, DIMENSION(:), POINTER :: p_sconfiguration
!</output>

!</subroutine>

  INTEGER :: ilength,isection,ivalue,ientry,icount
  CHARACTER, DIMENSION(:), POINTER :: p_sbuf
  
    IF (rparlist%isectionCount .EQ. 0) THEN
      PRINT *,'parlst_getStringRepresentation: Parameter list not initialised!'
      CALL sys_halt()
    END IF
  
    NULLIFY(p_sbuf)
    
    ! Number of characters in the buffer
    ilength = 0
    
    ! Loop through all sections
    DO isection = 1,rparlist%isectionCount
    
      ! Append the section name. May be empty for the unnamed section,
      ! which is always the first one.
      IF (isection .GT. 1) THEN
        ! Empty line before
        IF (ilength .GT. 0) CALL appendString(p_sbuf,ilength,'')
        CALL appendString(p_sbuf,ilength,&
          '['//TRIM(rparlist%p_Rsections(isection)%ssectionName)//']')
      END IF
        
      ! Loop through the values in the section
      DO ivalue = 1,rparlist%p_Rsections(isection)%iparamCount
        
        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        IF (icount .EQ. 0) THEN
          ! Write "name=value"
          CALL appendString(p_sbuf,ilength,&
            TRIM(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"="// &
            TRIM(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%sentry))
        ELSE
          ! Write "name(icount)="
          CALL appendString(p_sbuf,ilength,&
            TRIM(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//TRIM(sys_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          DO ientry = 1,icount
            CALL appendString(p_sbuf,ilength,&
              TRIM(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
                   p_Sentry(ientry)))
          END DO
        END IF
      
      END DO ! ivalue
    
    END DO ! isection
    
    ! Allocate a new character array with the correct size, copy p_sbuf to
    ! that ald release the old p_sbuf.
    ! Return NULL() if there is no data.
    NULLIFY(p_sconfiguration)
    IF (ilength .GT. 0) THEN
      ALLOCATE(p_sconfiguration(ilength))
      p_sconfiguration = p_sbuf(1:ilength)
    END IF
    
    ! Release our temp buffer
    IF (ASSOCIATED(p_sbuf)) DEALLOCATE(p_sbuf)

  CONTAINS
  
    ! Makes sure, the character buffer points to a character memory block of
    ! size nsize. If not, the block is reallocated to have that size.
    SUBROUTINE assumeBufSize(p_sconfig,nsize)
    
    CHARACTER, DIMENSION(:), POINTER :: p_sconfig
    INTEGER, INTENT(IN) :: nsize
    
    CHARACTER, DIMENSION(:), POINTER :: p_sconfignew
    
      IF (.NOT. ASSOCIATED(p_sconfig)) THEN
        ALLOCATE(p_sconfig(nsize))
      ELSE IF (SIZE(p_sconfig) .LT. nsize) THEN
        ALLOCATE(p_sconfignew(nsize))
        p_sconfignew(1:SIZE(p_sconfig)) = p_sconfig
        DEALLOCATE(p_sconfig)
        p_sconfig => p_sconfignew
      END IF
    
    END SUBROUTINE
    
    ! Appends sstring to the buffer p_sconfig, followed by a NEWLINE
    ! character. Reallocates memory if necessary.
    SUBROUTINE appendString(p_sconfig,iconfigLength,sstring)
    
    ! Pointer to character data
    CHARACTER, DIMENSION(:), POINTER :: p_sconfig
    
    ! In: Current length of data stream in p_sconfig.
    ! Out: New length of data stream in p_sconfig
    INTEGER, INTENT(INOUT) :: iconfigLength
    
    ! The string to be added.
    CHARACTER(LEN=*), INTENT(IN) :: sstring
    
      INTEGER :: nblocks,nblocksneeded,i
    
      ! How many memory blocks do we need for the current configuration?
      ! We work block-wise to prevent too often reallocation.
      IF (.NOT. ASSOCIATED(p_sconfig)) THEN
        nblocks = 0
      ELSE
        nblocks = SIZE(p_sconfig) / SYS_STRLEN
      END IF
      nblocksneeded = 1 + (iconfigLength+LEN(sstring)+1) / SYS_STRLEN
      IF (nblocksneeded .GT. nblocks) THEN
        CALL assumeBufSize(p_sconfig,nblocksneeded*SYS_STRLEN)
      END IF
      
      ! Append the data
      DO i=1,LEN(sstring)
        iconfigLength = iconfigLength+1
        p_sconfig(iconfigLength) = sstring(i:i)
      END DO
      
      ! Append NEWLINE as line-end character
      iconfigLength = iconfigLength+1
      p_sconfig(iconfigLength) = NEWLINE
    
    END SUBROUTINE
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE parlst_info (rparlist)
  
!<description>
  ! Prints the parameter list rparlist to the terminal.
!</description>
  
!<input> 
  ! The parameter list which is to be printed to the terminal.
  TYPE(t_parlist), INTENT(IN) :: rparlist
!</input>

!</subroutine>


  INTEGER :: isection,ivalue,ientry,icount
  
    IF (rparlist%isectionCount .EQ. 0) THEN
      PRINT *,'parlst_info: Parameter list not initialised!'
      CALL sys_halt()
    END IF
  
    ! Loop through all sections
    DO isection = 1,rparlist%isectionCount
    
      ! Append the section name. May be empty for the unnamed section,
      ! which is always the first one.
      IF (isection .GT. 1) THEN
        ! Empty line before
        CALL output_lbrk()
        CALL output_line('['//TRIM(rparlist%p_Rsections(isection)%ssectionName)//']')
      END IF
        
      ! Loop through the values in the section
      DO ivalue = 1,rparlist%p_Rsections(isection)%iparamCount
        
        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        IF (icount .EQ. 0) THEN
          ! Write "name=value"
          CALL output_line(&
            TRIM(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"="// &
            TRIM(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%sentry))
        ELSE
          ! Write "name(icount)="
          CALL output_line(&
            TRIM(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//TRIM(sys_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          DO ientry = 1,icount
            CALL output_line(&
              TRIM(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
                   p_Sentry(ientry)))
          END DO
        END IF
      
      END DO ! ivalue
    
    END DO ! isection
    
  END SUBROUTINE
  
END MODULE
