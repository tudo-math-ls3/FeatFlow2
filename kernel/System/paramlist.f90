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
!# 
!# Empty lines are ignored. When there is a '#' character in
!# a line not enclosed by apostrophes, the rest of the line is
!# ignored as a comment.
!#
!# The following routines can be used to maintain a parameter
!# list:
!# 1.) parlst_init 
!#      -> Initialises an empty parameter list
!#
!# 2.) parlst_readfromfile
!#     -> Reads the content of a .INI file into a parameter list.
!#
!# 3.) parlst_done
!#     -> Cleans up a parameter list, releases all allocated memory
!#        from the heap
!#
!# 4.) parlst_querysection
!#     -> Determines whether or not a section exists
!#
!# 5.) parlst_addsection
!#     -> Adds a new section
!#
!# 6.) parlst_queryvalue
!#     -> Determines whether or not a parameter exists
!#
!# 7.) parlst_getvalue_string
!#     -> Get the string value of a parameter from the parameter list
!#
!# 8.) parlst_addvalue
!#     -> Adds a new parameter to the parameter list
!#
!# 9.) parlst_setvalue
!#     -> Modifies the value of a parameter in the list
!# 
!# </purpose>
!##############################################################################

MODULE paramlist

  USE fsystem
  USE io
  
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
    
    ! A list of parameter valses. Each name contains PARLST_MLDATA
    ! characters. The length is identical to p_sparameters.
    CHARACTER(LEN=PARLST_MLDATA), DIMENSION(:), POINTER :: p_Svalues => NULL()
    
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

CONTAINS
  
  ! ***************************************************************************

  ! Internal subroutine: Initialise a newly created parameter section.
  
  SUBROUTINE parlst_initsection (rparlstSection,sname)
  
  TYPE(t_parlstSection), INTENT(INOUT) :: rparlstSection
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! Simply allocate the pointers with an empty list
  ALLOCATE(rparlstSection%p_Sparameters(PARLST_NPARSPERBLOCK))
  ALLOCATE(rparlstSection%p_Svalues(PARLST_NPARSPERBLOCK))
  
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
  
  INTEGER :: sz

  ! Pointers to new lists for replacing the old.
  CHARACTER(LEN=PARLST_MLNAME), DIMENSION(:), POINTER :: p_Sparameters 
  CHARACTER(LEN=PARLST_MLDATA), DIMENSION(:), POINTER :: p_Svalues 
  
  sz = MAX(SIZE(rparlstSection%p_Sparameters),inewsize)
  
  IF (SIZE(rparlstSection%p_Sparameters) .EQ. sz) RETURN ! nothing to do
  
  ! Allocate the pointers for the new lists
  ALLOCATE(p_Sparameters(sz))
  ALLOCATE(p_Svalues(sz))
  
  ! Copy the content of the old ones
  p_Sparameters(1:sz) = rparlstSection%p_Sparameters (1:sz)
  p_Svalues(1:sz) = rparlstSection%p_Svalues (1:sz)
  
  ! Throw away the old arrays, replace by the new ones
  DEALLOCATE(rparlstSection%p_Svalues)
  DEALLOCATE(rparlstSection%p_Sparameters)
  
  rparlstSection%p_Sparameters => p_Sparameters
  rparlstSection%p_Svalues => p_Svalues
  
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
    
    DEALLOCATE(rparlist%p_Rsections(i)%p_Svalues)
    DEALLOCATE(rparlist%p_Rsections(i)%p_Sparameters)
    
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
    STOP
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
  
  !</description>
  
!</function>

  ! local variables
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  exists = 0
  
  IF (sparameter .EQ. '') THEN
    PRINT *,'Empty parameter name!'
    STOP
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
  
  !</description>
  
!</function>

  ! local variables
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  exists = 0
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    STOP
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
  
!<subroutine>
  SUBROUTINE parlst_getvalue_string_indir (rsection, &
                            sparameter, svalue, sdefault)
  !<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, sdefault is returned.
  ! If sdefault is not given, an error will be thrown.
  
  !</description>
  
  !<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The parameter name.
  CHARACTER(LEN=*), INTENT(IN) :: sparameter

  ! Optional: A default value
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sdefault
  
  !</input>
  
  !<output>

  ! The value of the parametzer
  CHARACTER(LEN=*), INTENT(OUT) :: svalue
  
  !</output>

!</subroutine>

  ! local variables
  INTEGER :: i
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  IF (sparameter .EQ. '') THEN
    PRINT *,'Empty parameter name!'
    STOP
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
      PRINT *,'Parameter ',TRIM(paramname),'does not exist!'
      STOP
    END IF
  ELSE
    svalue = rsection%p_Svalues (i)
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_string_fetch (rsection, &
                            iparameter, svalue, iexists)
  !<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  ! If the value does not exist, sdefault is returned.
  ! If iexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If iexists is given, it will be set to YES if the parameter number
  ! iparameter exists, otherwise it will be set to NO and svalue=''.
  
  !</description>
  
  !<input>
    
  ! The section where to search for the parameter
  TYPE(t_parlstSection), INTENT(IN) :: rsection

  ! The number of the parameter.
  INTEGER, INTENT(IN) :: iparameter

  !</input>
  
  !<output>

  ! The value of the parameter
  CHARACTER(LEN=*), INTENT(OUT) :: svalue
  
  ! Optional: Parameter existance check
  ! Is set to YES/NO, depending on whether the parameter exists.
  INTEGER, INTENT(OUT), OPTIONAL :: iexists
  
  !</output>

!</subroutine>

  ! Check if iparameter is out of bounds. If yes, probably
  ! throw an error.
  
  IF ((iparameter .LT. 0) .OR. (iparameter .GT. rsection%iparamCount)) THEN
  
    IF (.NOT. PRESENT(iexists)) THEN 
      PRINT *,'Error. Parameter ',iparameter,' does not exist!'
      STOP
    ELSE
      svalue = ''
      iexists = NO
      RETURN
    END IF
  
  END IF
  
  ! Get the parameter value.
  svalue = rsection%p_Svalues (iparameter)
  IF (PRESENT(iexists)) iexists = YES

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_getvalue_string_direct (rparlist, ssectionName, &
                            sparameter, svalue, sdefault)
  !<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, sdefault is returned.
  ! If sdefault is not given, an error will be thrown.
  
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
    STOP
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Section not found'
    RETURN
  END IF

  ! Get the value
  IF (PRESENT(sdefault)) THEN
    CALL parlst_getvalue_string_indir (p_rsection, sparameter, svalue, sdefault)
  ELSE
    CALL parlst_getvalue_string_indir (p_rsection, sparameter, svalue)
  END IF  

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE parlst_addvalue_indir (rsection, sparameter, svalue)
  
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
  rsection%p_Svalues(rsection%iparamCount) = svalue

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_addvalue_direct (rparlist, ssectionName, sparameter, svalue)
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
  
  !</input>

!</subroutine>

  ! local variables
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    STOP
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Section not found'
    RETURN
  END IF

  ! Add the parameter 
  
  CALL parlst_addvalue_indir (p_rsection, sparameter, svalue)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE parlst_setvalue_fetch (rsection, iparameter, svalue, iexists)
  
  !<description>
  
  ! Modifies the value of a parameter in the section rsection.
  ! The value of parameter iparameter in the section rsection is modified.
  ! If iexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If iexists is given, it will be set to YES if the parameter number
  ! iparameter exists and was modified, otherwise it will be set to NO.
  
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
  
  !</input>

  !<output>

  ! Optional parameter. Is set to YES/NO, depending on whether
  ! the parameter exists.
  INTEGER, INTENT(OUT), OPTIONAL :: iexists

  !</output>

!</subroutine>

  ! Check if iparameter is out of bounds. If yes, probably
  ! throw an error.
  
  IF ((iparameter .LT. 0) .OR. (iparameter .GT. rsection%iparamCount)) THEN
  
    IF (.NOT. PRESENT(iexists)) THEN 
      PRINT *,'Error. Parameter ',iparameter,' does not exist!'
      STOP
    ELSE
      iexists = NO
      RETURN
    END IF
  
  END IF

  rsection%p_Svalues (iparameter) = svalue
  IF (PRESENT(iexists)) iexists = YES

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE parlst_setvalue_indir (rsection, sparameter, svalue)
  
  !<description>
  
  ! Modifies the value of a parameter in the section rsection.
  ! If the parameter does not exist, an error is thrown.
  
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
  
  !</input>

!</subroutine>

  ! local variables
  INTEGER :: i
  CHARACTER(LEN=PARLST_MLNAME) :: paramname
  
  ! Create the upper-case parameter name
  paramname = ADJUSTL(sparameter)
  CALL parlst_toupper (paramname)

  ! Get the parameter position
  i = parlst_queryvalue_indir (rsection, paramname)
  
  IF (i .EQ. 0) THEN
    PRINT *,'Parameter ',paramname,' does not exist, cannot be modified!'
    STOP
  ELSE 
    rsection%p_Svalues (i) = svalue
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE parlst_setvalue_direct (rparlist, ssectionName, sparameter, svalue)
  !<description>
  
  ! Modifies the value of a parameter in the section with name ssectionName
  ! in the parameter list rparlist.
  ! If the parameter does not exist, it's created.
  
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
  
  !</input>

!</subroutine>

  ! local variables
  TYPE(t_parlstSection), POINTER :: p_rsection
  
  ! Cancel if the list is not initialised.
  IF (rparlist%isectionCount .EQ. 0) THEN
    PRINT *,'Parameter list not initialised!'
    STOP
  END IF
  
  ! Get the section
  CALL parlst_querysection(rparlist, ssectionName, p_rsection) 
  IF (.NOT. ASSOCIATED(p_rsection)) THEN
    PRINT *,'Section not found'
    RETURN
  END IF

  ! Set the parameter 
  
  CALL parlst_setvalue_indir (p_rsection, sparameter, svalue)

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
  
  SUBROUTINE parlst_parseline (sdata, ityp, ilinenum, ssecname, sparamname, svalue)
  
  ! The line to be parsed
  CHARACTER(LEN=*), INTENT(IN) :: sdata
  
  ! The typ of the line
  INTEGER, INTENT(OUT) :: ityp

  ! Line number
  INTEGER, INTENT(IN) :: ilinenum
  
  ! Section name, if it's a section
  CHARACTER(LEN=*), INTENT(INOUT) :: ssecname
  
  ! Parameter name, if it's a parameter
  CHARACTER(LEN=*), INTENT(INOUT) :: sparamname
  
  ! Parameter value, if it's a parameter
  CHARACTER(LEN=*), INTENT(INOUT) :: svalue
  
  ! local variables
  INTEGER :: i,ltr
  CHARACTER(LEN=PARLST_LENLINEBUF) :: sbuf
  
    ityp = 0
    
    ! Do we have data in sdata?
    IF (sdata .EQ. '') RETURN
    
    ! Copy the input string - left adjusted - and get the string length
    sbuf = ADJUSTL(sdata)
    
    ! Do we start with '[' and end with ']'?
    IF (sbuf(1:1) .EQ. "[") THEN
    
      ! Find the final ']'.
      DO ltr = 1,LEN(sbuf)
        IF (sbuf(ltr:ltr) .EQ. "]") EXIT
      END DO
      
      IF (sbuf(ltr:ltr) .NE. ']') THEN
        PRINT *,'Wrong syntax of section name. Line ',ilinenum,':'
        PRINT *,sbuf
        STOP
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
      
      ! Is there a '=' sign?
      i = INDEX(sbuf(1:ltr),'=')
      IF (i .EQ. 0) THEN
        PRINT *,'Invalid parameter syntax. Line ',ilinenum,':'
        STOP
      END IF
    
      ityp = 2
      
      ! Get the name of the parameter
      sparamname = ADJUSTL(sbuf(1:i-1))
      
      ! Get the parameter value
      svalue = ADJUSTL(sbuf(i+1:ltr))
    
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
  INTEGER :: iunit,ios,isbuflen,ityp,ilinenum
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
    STOP
  END IF
  
  ! Start adding parameters to the unnamed section
  p_currentsection => rparlist%p_Rsections(1)
  
  ! Read all lines from the file
  ios = 0
  ilinenum = 0
  DO WHILE (ios .EQ. 0) 
    
    ! Read a line from the file into sbuf
    CALL parlst_readlinefromfile (iunit, sdata, isbuflen, ios)
    ilinenum = ilinenum + 1
    
    IF (isbuflen .NE. 0) THEN
    
      ! Parse the line
      CALL parlst_parseline (sdata, ityp, ilinenum, ssectionname, sparname, svalue)  
      
      SELECT CASE (ityp)
      CASE (1)
        ! A new section name. Add a section, set the current section
        ! to the new one.
        CALL parlst_addsection (rparlist, ssectionname)
        p_currentsection => rparlist%p_Rsections(rparlist%isectionCount)
        
      CASE (2)
        ! A new parameter. Add it to the current section.
        CALL parlst_addvalue (p_currentsection, sparname, svalue)
        
      ! Other cases: comment.
      END SELECT
    
    END IF
  
  END DO
  
  ! Close the file, finish.
  CLOSE (iunit)

  END SUBROUTINE
  
END MODULE
