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
!# The subvariables feature
!# ========================
!# The parameter list alwo allows to specify variables as subvariables
!# of other variables. Take a look at the following example:
!# 
!# -------------------------snip------------------------------
!# NLMIN=1
!# ACOMPLEXXTRING=! %{NLMIN} ! %{imainelement} ! %{SECTION1.ilist:3}
!# imainelement=%{SECTION1.ielement}
!# 
!# [SECTION1]
!# ielement=5
!# 
!# ilist(4)=
!#   abc
!#   def
!#   ghi
!#   jkl
!# -------------------------snip------------------------------
!# When reading the file, the parameter "ACOMPLEXXTRING" is automatically
!# expanded to the value "! 1 ! 5 ! ghi" by using other variables from
!# the parameter file. The following syntax is allowed here to refer
!# to other parameters:
!#
!#  %{NAME}               - A variable from the unnamed section
!#  %{NAME:idx}           - Value number idx of variable NAME 
!#                          from the unnamed section
!#  %{NAME:idx}           - Variable NAME from section SECTION
!#  %{SECTION.NAME:idx}   - Value number idx of variable NAME 
!#                          from section SECTION
!#
!# The environment variable feature
!# ================================
!# DAT files may refer to environment variables. Example:
!# -------------------------snip------------------------------
!# NLMIN=1
!# NLMAX=$NLMAXENV
!# -------------------------snip------------------------------
!# In this case, it's assumed that $NLMAXENV is an environment
!# variable of the system where the program runs. So when NLMAX
!# is requested, the value of the environment variable $NLMAXENV
!# is returned.
!#
!# The subfile feature
!# ===================
!# An INI file may contain references to subfiles. Subfiles must
!# be specified at the beginning of an INI file with the following 
!# syntax:
!# -------------------------snip------------------------------
!# # this is a data file which imports a couple of subfiles
!#
!# simportdatafiles(4) = 
!#   "subfilename1.ini"
!#   "subfilename2.ini"
!#   "subfilename3.ini"
!#   "subfilename4.ini"
!#
!# # The rest of the file overwrites the data in the subfiles.
!#
!# parameter1 = data1
!# parameter2 = data2
!#
!# # now a section follows; if one of the child data files above
!# # also contains that section and the parameter 'parameter3',
!# # that data is overwritten here!
!#
!# [Section1]
!# parameter3 = 'this string replaces parameter3 from the child-files'
!# ...
!# -------------------------snip------------------------------
!#
!# parlst_readfromfile will at first read all subfiles and then
!# evaluate the data in the main file. Data in the main file will
!# overwrite data in the subfiles. Subfiles may contain references
!# to other subfiles. That way, a 'master.ini' file may define
!# a couple of files which are read at first and then overwrite
!# some parameters.
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
!#      parlst_getvalue_single
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
!# Auxiliary routines
!#
!# 1.) parlst_readfromsinglefile
!#     -> Reads a single INI file without checking for subfiles.
!#
!# 2.) parlst_expandEnvVariables
!#     -> expands all environment variables in the given
!#        string sbuffer.
!#
!# 3.) parlst_expandsubvars
!#     -> expands all subvariables in the parameter list
!#
!# </purpose>
!##############################################################################

module paramlist

  use fsystem
  use io
  use genoutput
  
  implicit none

!<constants>

  !<constantblock>
  
  ! Maximum length of a section name.
  integer, parameter :: PARLST_MLSECTION = 64

  ! Maximum length of parameter names: 32 characters
  integer, parameter :: PARLST_MLNAME = 32

  ! Maximum length of parameter data: 256 characters
  integer, parameter :: PARLST_MLDATA = 256
  
  ! Minimum number of free parameter 'slots' per parameter section.
  ! If there are too many parameters in a parameter section, the
  ! structure is dynamically extended in terms of PARLST_NPARSPERBLOCK
  ! entries.
  integer, parameter :: PARLST_NPARSPERBLOCK = 32

  ! Minimum number of parameter sections.
  ! If there are too many parameter sections in a parameter block, the
  ! structure is dynamically extended in terms of PARLST_NSECTIONS
  ! entries.
  integer, parameter :: PARLST_NSECTIONS = 8

  ! Maximum length of a line in a INI file. Lines longer than this
  ! are truncated.
  integer, parameter :: PARLST_LENLINEBUF = 1024
  
  ! Comment character
  character, parameter :: PARLST_COMMENT = "#"

  !</constantblock>

!</constants>

!<types>
  
  !<typeblock>
  
  ! This structure realises a value associated to a parameter name.
  ! A value consists of one string or an array of strings.
  
  type t_parlstValue
  
    ! Number of strings. If set to 0, the value consists of one
    ! string, to be found in svalue. If > 0, there are nsize
    ! strings to be found in p_Sentry.
    integer :: nsize = 0
    
    ! Single string; contains the value in case nsize=0
    character(LEN=PARLST_MLDATA) :: sentry = ''
    
    ! Array of strings in case nsize>0
    character(LEN=PARLST_MLDATA), dimension(:), pointer :: p_Sentry => null()
  
  end type
  
  !</typeblock>
  
  !<typeblock>
  
  ! This structure realises a parameter section. It contains an
  ! array with parameter names and an array with parameter values
  ! to these names. The arrays are dynamically allocated. 
  
  type t_parlstSection
  
    ! The name of the section.
    character(LEN=PARLST_MLSECTION) :: ssectionName = ''
    
    ! Actual number of parameters in this section.
    integer :: iparamCount = 0
    
    ! A list of parameter names. Each name contains PARLST_MLNAME
    ! characters.
    character(LEN=PARLST_MLNAME), dimension(:), pointer :: p_Sparameters => null()
    
    ! A list of t_parlstValue structures corresponding to the parameters
    ! in p_Sparameters.
    type(t_parlstValue), dimension(:), pointer :: p_Rvalues
    
  end type
  
  !</typeblock>
  
  !<typeblock>
  
  ! This structure realises a parameter list. Parameters can be read into
  ! it from a file. Parameters can be obtained from the structure using
  ! the query/get routines.
  
  type t_parlist
  
    ! Actual number of sections in the parameter list. There's at least
    ! one section - the unnamed section. If this value is =0, the parameter
    ! list is not initialised.
    integer :: isectionCount = 0
    
    ! A list of sections. The first section is always the unnamed section.
    type(t_parlstSection), dimension(:), pointer :: p_Rsections => null()
    
  end type
  
  !</typeblock>

!</types>

  private :: parlst_initsection, parlst_reallocsection, parlst_realloclist
  private :: parlst_fetchparameter,parlst_readlinefromfile,parlst_parseline
  
  interface parlst_queryvalue
    module procedure parlst_queryvalue_direct
    module procedure parlst_queryvalue_indir
  end interface

  interface parlst_querysubstrings
    module procedure parlst_querysubstrings_direct
    module procedure parlst_querysubstrings_indir
  end interface

  interface parlst_addvalue
    module procedure parlst_addvalue_direct
    module procedure parlst_addvalue_indir
  end interface

  interface parlst_setvalue
    module procedure parlst_setvalue_fetch
    module procedure parlst_setvalue_indir
    module procedure parlst_setvalue_direct
  end interface

  interface parlst_getvalue_string
    module procedure parlst_getvalue_string_fetch
    module procedure parlst_getvalue_string_indir
    module procedure parlst_getvalue_string_direct
  end interface

  interface parlst_getvalue_int
    module procedure parlst_getvalue_int_fetch
    module procedure parlst_getvalue_int_indir
    module procedure parlst_getvalue_int_direct
  end interface

  interface parlst_getvalue_single
    module procedure parlst_getvalue_single_fetch
    module procedure parlst_getvalue_single_indir
    module procedure parlst_getvalue_single_direct
  end interface

  interface parlst_getvalue_double
    module procedure parlst_getvalue_double_fetch
    module procedure parlst_getvalue_double_indir
    module procedure parlst_getvalue_double_direct
  end interface

contains
  
  ! ***************************************************************************

  ! Internal subroutine: Initialise a newly created parameter section.
  
  subroutine parlst_initsection (rparlstSection,sname)
  
  type(t_parlstSection), intent(INOUT) :: rparlstSection
  character(LEN=*), intent(IN) :: sname
  
  ! Simply allocate the pointers with an empty list
  allocate(rparlstSection%p_Sparameters(PARLST_NPARSPERBLOCK))
  allocate(rparlstSection%p_Rvalues(PARLST_NPARSPERBLOCK))
  
  ! and set the section name
  rparlstSection%ssectionName = sname
  
  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Reallocate a section.
  ! This increases the size of a parameter section by reallocation of the
  ! arrays.
  
  subroutine parlst_reallocsection (rparlstSection, inewsize)
  
  ! The section to reallocate.
  type(t_parlstSection), intent(INOUT) :: rparlstSection
  
  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the section should be able to handle.
  integer, intent(IN) :: inewsize
  
  ! local variables
  
  integer :: sz,oldsize

  ! Pointers to new lists for replacing the old.
  character(LEN=PARLST_MLNAME), dimension(:), pointer :: p_Sparameters 
  type(t_parlstValue), dimension(:), pointer :: p_Rvalues
  
  oldsize = size(rparlstSection%p_Sparameters)
  sz = max(oldsize,inewsize)
  
  if (size(rparlstSection%p_Sparameters) .eq. sz) return ! nothing to do
  
  ! Allocate the pointers for the new lists
  allocate(p_Sparameters(sz))
  allocate(p_Rvalues(sz))
  
  ! Copy the content of the old ones
  p_Sparameters(1:oldsize) = rparlstSection%p_Sparameters (1:oldsize)
  p_Rvalues(1:oldsize) = rparlstSection%p_Rvalues (1:oldsize)
  
  ! Throw away the old arrays, replace by the new ones
  deallocate(rparlstSection%p_Rvalues)
  deallocate(rparlstSection%p_Sparameters)
  
  rparlstSection%p_Sparameters => p_Sparameters
  rparlstSection%p_Rvalues => p_Rvalues
  
  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Release a section.
  ! Removes all temporary memory that is allocated by a section.
  
  subroutine parlst_releasesection (rparlstSection)
  
  ! The section to release.
  type(t_parlstSection), intent(INOUT) :: rparlstSection
  
  ! local variables
  integer :: i
  
  ! Loop through all values in the current section if there is
  ! an array-value. Release them.
  do i=size(rparlstSection%p_Rvalues),1,-1
    if (rparlstSection%p_Rvalues(i)%nsize .gt. 0) then
      deallocate(rparlstSection%p_Rvalues(i)%p_Sentry)
    end if
  end do
  
  ! Remove the content of the section.
  deallocate(rparlstSection%p_Rvalues)
  deallocate(rparlstSection%p_Sparameters)
  rparlstSection%iparamCount = 0
  
  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Reallocate the section list
  ! This increases the size of a section list by reallocation of the
  ! arrays.
  
  subroutine parlst_realloclist (rparlist, inewsize)
  
  ! The section list to reallocate.
  type(t_parlist), intent(INOUT) :: rparlist
  
  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the section should be able to handle.
  integer, intent(IN) :: inewsize
  
  ! local variables
  
  integer :: sz

  ! Pointers to new lists for replacing the old.
  type(t_parlstSection), dimension(:), pointer :: p_Rsections
  
  ! Allocate the pointers for the new lists
  allocate(p_Rsections(inewsize))
  
  sz = min(size(rparlist%p_Rsections),inewsize)
  
  ! Copy the content of the old ones
  p_Rsections(1:sz) = rparlist%p_Rsections (1:sz)
  
  ! Throw away the old arrays, replace by the new ones
  deallocate(rparlist%p_Rsections)
  
  rparlist%p_Rsections => p_Rsections
  
  end subroutine
  
  ! ***************************************************************************
  
  ! Internal subroutine: Convert a character to upper case.
  
  subroutine parlst_toupper (str) 
  
  ! The string that is to make uppercase
  character(LEN=*), intent(INOUT) :: str
  
!  CHARACTER(LEN=26), PARAMETER :: lowc = 'abcdefghijklmnopqrstuvwxyz'
!  CHARACTER(LEN=26), PARAMETER :: upc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!  
!  INTEGER :: i
!  
!  i = INDEX(lowc,c)
!  IF (i .NE. 0) THEN
!    c = upc(i:i)
!  END IF

  integer, parameter :: up2low = iachar("a") - iachar("A")
  integer :: i
  character    :: c
  
  do i=1,len(str)
    c = str(i:i)
    if ((c .ge. "a") .and. (c .le. "z")) then
      str(i:i) = achar (iachar(c) - up2low)
    end if
  end do
  
  end subroutine
  
  ! ***************************************************************************
  
  ! Internal subroutine: Search in a section for a parameter
  ! and return the index - or 0 if the parameter does not exist.

  subroutine parlst_fetchparameter(rsection, sname, iparamnum) 

  ! The section.
  type(t_parlstSection), intent(IN) :: rsection
  
  ! The parameter name to look for. Must be uppercase.
  character(LEN=*), intent(IN) :: sname
  
  ! The number of the parameter in the list or 0 if it does not exist.
  integer, intent(OUT) :: iparamnum
  
  ! local variables
  integer :: i
  
  iparamnum = 0
  
  ! If the parameter list is empty, the section does not exist for sure
  if (rsection%iparamCount .eq. 0) return
  
  ! Loop through all sections to see if the section exists
  do i=1,rsection%iparamCount
    if (rsection%p_Sparameters(i) .eq. sname) then
      iparamnum = i
      return
    end if
  end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine parlst_init (rparlist)
  
!<description>
  
  ! This routine initialises a parameter list. It must be applied to a
  ! parameter list structure before doing anything to it, just to initialise.
  
!</description>
  
!<inputoutput>
  
  ! The parameter list to initialise.
  type(t_parlist), intent(INOUT) :: rparlist
  
!</inputoutput>
  
!</subroutine>

  ! Set the section-count to 1.
  rparlist%isectionCount = 1
  
  ! Allocate a first set of sections
  allocate(rparlist%p_Rsections(PARLST_NSECTIONS))
  
  ! Initialise the first section - it's the unnamed one.
  call parlst_initsection (rparlist%p_Rsections(1),'')
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine parlst_clear (rparlist)
  
!<description>
  ! This routine cleans up a parameter list. All parameters in rparlist are
  ! removed.
!</description>
  
!<inputoutput>
  ! The parameter list to clean up.
  type(t_parlist), intent(INOUT) :: rparlist
!</inputoutput>
  
!</subroutine>

    ! Clean up = done+reinit. We make that simple here...
    call parlst_done (rparlist)
    call parlst_init (rparlist)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine parlst_done (rparlist)
  
!<description>
  
  ! This routine releases a parameter list. All memory allocated by the
  ! parameter list is released.
  
!</description>
  
!<inputoutput>
  
  ! The parameter list to release.
  type(t_parlist), intent(INOUT) :: rparlist
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

  ! Probably nothing to do
  if (rparlist%isectionCount .eq. 0) return

  ! Loop through the parameter lists and release the content
  do i=rparlist%isectionCount,1,-1
    call parlst_releasesection (rparlist%p_Rsections(i))
  end do

  ! Release all sections
  deallocate(rparlist%p_Rsections)
  
  ! Mark the structure as 'empty', finish
  rparlist%isectionCount = 0

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine parlst_querysection(rparlist, sname, p_rsection) 

!<description>

  ! Searches for a section and return a pointer to it -
  ! or NULL() of the section does not exist.
  
!</description>

!<input>

  ! The parameter list to scan for the section.
  type(t_parlist), intent(IN) :: rparlist
  
  ! The section name to look for. 
  character(LEN=*), intent(IN) :: sname
  
!</input>
  
!<output>
  
  ! A pointer to the section.
  type(t_parlstSection), pointer :: p_rsection
  
!</output>
  
!</subroutine>
  
  ! local variables
  integer :: i
  character(LEN=PARLST_MLSECTION) :: sectionname
  
  nullify(p_rsection)
  
  ! If the parameter list is empty, the section does not exist for sure
  if (rparlist%isectionCount .eq. 0) return
  
  ! If the section name is '', return a pointer to the first section.
  if (sname .eq. '') then
    p_rsection => rparlist%p_Rsections(1)
    return
  end if
  
  ! Create the upper-case section name
  sectionname = adjustl(sname)
  call parlst_toupper (sectionname)

  ! Loop through all sections to see if the section exists
  do i=1,rparlist%isectionCount
    if (rparlist%p_Rsections(i)%ssectionName .eq. sectionname) then
      p_rsection => rparlist%p_Rsections(i)
      return
    end if
  end do

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine parlst_addsection (rparlist, sname)
  
!<description>
  
  ! Adds a section with the name sname to the list of sections in the
  ! parameter list rparlist. The name must NOT contain brackets ('[',']')
  ! in front and at the end!
  
!</description>

!<inputoutput>
  
  ! The parameter list where to add the section.
  type(t_parlist), intent(INOUT) :: rparlist
  
!</inputoutput>

!<input>
  
  ! The section name to add - without brackets in front and at the end!
  character(LEN=*), intent(IN) :: sname
  
!</input>
  
!</subroutine>

  ! local variables
  character(LEN=PARLST_MLSECTION) :: sectionname
  
  ! Cancel if the list is not initialised.
  if (rparlist%isectionCount .eq. 0) then
    call output_line ('Parameter list not initialised!', &
            OU_CLASS_ERROR,OU_MODE_STD,'parlst_addsection')
    call sys_halt()
  end if
  
  ! Create the upper-case section name
  sectionname = adjustl(sname)
  call parlst_toupper (sectionname)
  
  ! Add a new section - reallocate the section list if necessary
  if (rparlist%isectionCount .eq. size(rparlist%p_Rsections)) then
    call parlst_realloclist (rparlist, size(rparlist%p_Rsections)+PARLST_NSECTIONS)
  end if
  rparlist%isectionCount = rparlist%isectionCount + 1
  
  ! Initialise the new section.
  call parlst_initsection(rparlist%p_Rsections(rparlist%isectionCount),sectionname)

  end subroutine
  
  ! ***************************************************************************
  
!<function>

  integer function parlst_queryvalue_indir (rsection, sparameter) &
               result (exists)
          
!<description>
  ! Checks whether a parameter sparameter exists in the section rsection.
!</description>
  
!<result>
  ! The index of the parameter in the section ssection or =0, if the
  ! parameter does not exist within the section.
!</result>

!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The parameter name to search for.
  character(LEN=*), intent(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  character(LEN=PARLST_MLNAME) :: paramname
  
  exists = 0
  
  if (sparameter .eq. '') then
    call output_line ('Empty parameter name!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_queryvalue_indir')
    call sys_halt()
  end if
  
  ! Create the upper-case parameter name
  paramname = adjustl(sparameter)
  call parlst_toupper (paramname)
  
  ! Get the parameter index into 'exists', finish.
  call parlst_fetchparameter(rsection, paramname, exists)
  
  end function
  
  ! ***************************************************************************
  
!<function>

  integer function parlst_queryvalue_direct (rparlist, ssectionName, sparameter) &
               result (exists)
          
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
  type(t_parlist), intent(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name to search for.
  character(LEN=*), intent(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  type(t_parlstSection), pointer :: p_rsection
  
  exists = 0
  
  ! Cancel if the list is not initialised.
  if (rparlist%isectionCount .eq. 0) then
    call output_line ('Parameter list not initialised!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_queryvalue_direct')
    call sys_halt()
  end if
  
  ! Get the section
  call parlst_querysection(rparlist, ssectionName, p_rsection) 
  if (.not. associated(p_rsection)) then
    call output_line ('Section not found: '//trim(ssectionName), &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_queryvalue_direct')
    return
  end if
  
  ! Search for the parameter
  exists = parlst_queryvalue_indir (p_rsection, sparameter)

  end function

  ! ***************************************************************************
  
!<function>

  integer function parlst_querysubstrings_indir (rsection, sparameter) &
               result (iresult)
          
!<description>
  ! Returns the number of substrings of a parameter.
!</description>
  
!<result>
  ! The number of substrings of parameter sparameter in section rsection.
!</result>

!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The parameter name to search for.
  character(LEN=*), intent(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  integer :: idx
  character(LEN=PARLST_MLNAME) :: paramname
  
  if (sparameter .eq. '') then
    call output_line ('Empty parameter name!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_querysubstrings_indir')
    call sys_halt()
  end if
  
  ! Create the upper-case parameter name
  paramname = adjustl(sparameter)
  call parlst_toupper (paramname)
  
  ! Get the parameter index into 'idx', finish.
  call parlst_fetchparameter(rsection, paramname, idx)
  
  ! Return number of substrings
  if (idx .eq. 0) then
    iresult = 0
  else
    iresult = rsection%p_Rvalues(idx)%nsize
  end if
  
  end function
  
  ! ***************************************************************************
  
!<function>

  integer function parlst_querysubstrings_direct (rparlist, ssectionName, sparameter) &
               result (iresult)
          
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
  type(t_parlist), intent(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name to search for.
  character(LEN=*), intent(IN) :: sparameter
  
!</input>
  
!</function>

  ! local variables
  integer :: idx
  type(t_parlstSection), pointer :: p_rsection
  
  ! Cancel if the list is not initialised.
  if (rparlist%isectionCount .eq. 0) then
    call output_line ('Parameter list not initialised!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_querysubstrings_direct')
    call sys_halt()
  end if
  
  ! Get the section
  call parlst_querysection(rparlist, ssectionName, p_rsection) 
  if (.not. associated(p_rsection)) then
    call output_line ('Section not found: '//trim(ssectionName), &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_querysubstrings_direct')
    return
  end if
  
  ! Get the parameter index
  idx = parlst_queryvalue_indir (p_rsection, sparameter)

  ! Return number of substrings
  if (idx .eq. 0) then
    iresult = 0
  else
    iresult = p_rsection%p_Rvalues(idx)%nsize
  end if

  end function

  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_string_indir (rsection, sparameter, svalue, &
                                           sdefault, isubstring)
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
  ! When omitting isubstring, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  character(LEN=*), intent(IN), optional :: sdefault
  
  ! OPTIONAL: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  integer, intent(IN), optional :: isubstring
  
!</input>
  
!<output>

  ! The value of the parameter
  character(LEN=*), intent(OUT) :: svalue
  
!</output>

!</subroutine>

  ! local variables
  integer :: i,isub
  character(LEN=PARLST_MLNAME) :: paramname
  
  if (sparameter .eq. '') then
    call output_line ('Empty parameter name!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_getvalue_string_indir')
    call sys_halt()
  end if

  ! Create the upper-case parameter name
  paramname = adjustl(sparameter)
  call parlst_toupper (paramname)
  
  ! Get the parameter index into 'exists', finish.
  call parlst_fetchparameter(rsection, paramname, i)
  
  if (i .eq. 0) then
    if (present(sdefault)) then
      svalue = sdefault
    else
      call output_line ('Parameter not found: '//trim(paramname), &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_getvalue_string_indir')
      call sys_halt()
    end if
  else
    ! Depending on isubstring, return either the 'headline' or one
    ! of the substrings.
    isub = 0
    if (present(isubstring)) isub = isubstring
  
    if ((isub .le. 0) .or. (isub .gt. rsection%p_Rvalues(i)%nsize)) then
      svalue = rsection%p_Rvalues(i)%sentry
    else
      svalue = rsection%p_Rvalues(i)%p_Sentry(isub)
    end if
  end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_string_fetch (rsection, iparameter, svalue,&
                                           bexists, isubstring)
!<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  ! If bexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  !
  ! If bexists is given, it will be set to TRUE if the parameter number
  ! iparameter exists, otherwise it will be set to FALSE and svalue=''.
  !  
  ! If the value is an array of strings, the optional parameter isubstring>=0
  ! allows to specify the number of the substring to be returned; 
  ! isubstring=0 returns the value directly
  ! behind the '=' sign in the line of the parameter, isubstring>0 returns
  ! the array-entry in the lines below the parameter.
  !
  ! When omitting isubstring, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The number of the parameter.
  integer, intent(IN) :: iparameter

  ! OPTIONAL: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  integer, intent(IN), optional :: isubstring

!</input>
  
!<output>

  ! The value of the parameter
  character(LEN=*), intent(OUT) :: svalue
  
  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(OUT), optional :: bexists
  
!</output>

!</subroutine>

  integer :: isub

  ! Check if iparameter is out of bounds. If yes, probably
  ! throw an error.
  
  if ((iparameter .lt. 0) .or. (iparameter .gt. rsection%iparamCount)) then
  
    if (.not. present(bexists)) then 
      call output_line ('Error. Parameter '//trim(sys_siL(iparameter,10))//&
          ' does not exist!', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_getvalue_string_fetch')
      call sys_halt()
    else
      svalue = ''
      bexists = .false.
      return
    end if
  
  end if
  
  ! Get the parameter value.
  ! Depending on isubstring, return either the 'headline' or one
  ! of the substrings.
  isub = 0
  if (present(isubstring)) isub = isubstring

  if ((isub .le. 0) .or. &
      (isub .gt. rsection%p_Rvalues(iparameter)%nsize)) then
    svalue = rsection%p_Rvalues(iparameter)%sentry
  else
    svalue = rsection%p_Rvalues(iparameter)%p_Sentry(isub)
  end if
  
  if (present(bexists)) bexists = .true.

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_string_direct (rparlist, ssectionName, &
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
  ! When omitting isubstring, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The parameter list.
  type(t_parlist), intent(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  character(LEN=*), intent(IN), optional :: sdefault
  
  ! OPTIONAL: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  integer, intent(IN), optional :: isubstring

!</input>
  
!<output>

  ! The value of the parameter
  character(LEN=*), intent(OUT) :: svalue
  
!</output>

!</subroutine>

  ! local variables
  type(t_parlstSection), pointer :: p_rsection
  
  ! Cancel if the list is not initialised.
  if (rparlist%isectionCount .eq. 0) then
    call output_line ('Parameter list not initialised!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_getvalue_string_fetch')
    call sys_halt()
  end if
  
  ! Get the section
  call parlst_querysection(rparlist, ssectionName, p_rsection) 
  if (.not. associated(p_rsection)) then
    if (present(sdefault)) then
      svalue = sdefault
      return
    else
      call output_line ('Section not found: '//trim(ssectionName), &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_getvalue_string_fetch')
      call sys_halt()
    end if
  end if

  ! Get the value
  call parlst_getvalue_string_indir (p_rsection, sparameter, svalue, sdefault,&
                                     isubstring)

  ! Expand any environment variable
  svalue = parlst_expandEnvVariables(svalue)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_single_indir (rsection, sparameter, fvalue, &
                                           fdefault, iarrayindex)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, idefault is returned.
  ! If idefault is not given, an error will be thrown.
  !
  ! If the value is an array of singles, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the single to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  real(SP), intent(IN), optional :: fdefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  real(SP), intent(OUT) :: fvalue
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  if (present(fdefault)) then
    write (sdefault,'(E17.10E2)') fdefault
    call parlst_getvalue_string_indir (rsection, sparameter, svalue, &
                                       sdefault, iarrayindex)
  else
    call parlst_getvalue_string_indir (rsection, sparameter, svalue, &
                                       isubstring=iarrayindex)
  end if

  fvalue = sys_Str2Single(svalue,'(E17.10E2)')
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine parlst_getvalue_single_fetch (rsection, iparameter, fvalue, &
                                           bexists, iarrayindex)

!<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  !
  ! If bexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If bexists is given, it will be set to TRUE if the parameter number
  ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
  !
  ! If the value is an array of singles, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the single to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.

!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The number of the parameter.
  integer, intent(IN) :: iparameter
  
  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  real(SP), intent(OUT) :: fvalue
  
  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(OUT), optional :: bexists
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: svalue
  
  svalue = '0.0E0'
  call parlst_getvalue_string_fetch (rsection, iparameter, svalue, &
                                     bexists, iarrayindex)
  
  fvalue = sys_Str2Single(svalue,'(E17.10E2)')
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_single_direct (rparlist, ssectionName, &
                                            sparameter, fvalue, fdefault,&
                                            iarrayindex)
!<description>
  
  ! Returns the value of a parameter in the section ssection. If the
  ! value does not exist, ddefault is returned.  If fdefault is not
  ! given, an error will be thrown.
  !
  ! If the value is an array of singles, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the single to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.
!</description>
  
!<input>
    
  ! The parameter list.
  type(t_parlist), intent(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  real(SP), intent(IN), optional :: fdefault
  
  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  real(SP), intent(OUT) :: fvalue
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  if (present(fdefault)) then
    write (sdefault,'(E17.10E2)') fdefault
    call parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, sdefault, &
                                        iarrayindex)
  else
    call parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, &
                                        isubstring=iarrayindex)
  end if
  
  fvalue = sys_Str2Single(svalue,'(E17.10E2)')
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_double_indir (rsection, sparameter, dvalue, &
                                           ddefault, iarrayindex)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, idefault is returned.
  ! If idefault is not given, an error will be thrown.
  !
  ! If the value is an array of doubles, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the double to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  real(DP), intent(IN), optional :: ddefault
  
  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  real(DP), intent(OUT) :: dvalue
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  if (present(ddefault)) then
    write (sdefault,'(E27.19E3)') ddefault
    call parlst_getvalue_string_indir (rsection, sparameter, svalue, &
                                       sdefault, iarrayindex)
  else
    call parlst_getvalue_string_indir (rsection, sparameter, svalue, &
                                       isubstring=iarrayindex)
  end if

  dvalue = sys_Str2Double(svalue,'(E27.19E3)')
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine parlst_getvalue_double_fetch (rsection, iparameter, dvalue, &
                                           bexists, iarrayindex)

!<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  !
  ! If bexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  ! If bexists is given, it will be set to TRUE if the parameter number
  ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
  !
  ! If the value is an array of doubles, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the double to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.  
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The number of the parameter.
  integer, intent(IN) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  real(DP), intent(OUT) :: dvalue
  
  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(OUT), optional :: bexists
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: svalue
  
  svalue = '0.0E0'
  call parlst_getvalue_string_fetch (rsection, iparameter, svalue, &
                                     bexists, iarrayindex)
  
  dvalue = sys_Str2Double(svalue,'(E27.19E3)')
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_double_direct (rparlist, ssectionName, &
                                            sparameter, dvalue, ddefault,&
                                            iarrayindex)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, ddefault is returned.
  ! If ddefault is not given, an error will be thrown.
  !
  ! If the value is an array of doubles, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the double to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.

!</description>
  
!<input>
    
  ! The parameter list.
  type(t_parlist), intent(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  real(DP), intent(IN), optional :: ddefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  real(DP), intent(OUT) :: dvalue
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  if (present(ddefault)) then
    write (sdefault,'(E27.19E3)') ddefault
    call parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, sdefault, &
                                        iarrayindex)
  else
    call parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, &
                                        isubstring=iarrayindex)
  end if
  
  dvalue = sys_Str2Double(svalue,'(E27.19E3)')
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_int_indir (rsection, sparameter, ivalue, &
                                        idefault, iarrayindex)
!<description>
  
  ! Returns the value of a parameter in the section ssection.
  ! If the value does not exist, idefault is returned.
  ! If idefault is not given, an error will be thrown.
  !
  ! If the value is an array of integers, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the integer to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  integer, intent(IN), optional :: idefault
  
  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  integer, intent(OUT) :: ivalue
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  if (present(idefault)) then
    write (sdefault,*) idefault
    call parlst_getvalue_string_indir (rsection, sparameter, svalue, &
                                       sdefault, iarrayindex)
  else
    call parlst_getvalue_string_indir (rsection, sparameter, svalue, &
                                       isubstring=iarrayindex)
  end if
  
  read(svalue,*) ivalue

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_int_fetch (rsection, iparameter, ivalue, &
                                        bexists, iarrayindex)
!<description>
  
  ! Returns the value of a parameter in the section rsection.
  ! iparameter specifies the number of the parameter in section rsection.
  ! If bexists does not appear, an error is thrown if a nonexisting
  ! parameter is accessed.
  !
  ! If bexists is given, it will be set to TRUE if the parameter number
  ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
  !
  ! If the value is an array of integers, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the integer to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The section where to search for the parameter
  type(t_parlstSection), intent(IN) :: rsection

  ! The number of the parameter.
  integer, intent(IN) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  integer, intent(OUT) :: ivalue
  
  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(OUT), optional :: bexists
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: svalue
  
  svalue = '0'
  call parlst_getvalue_string_fetch (rsection, iparameter, svalue, &
                                     bexists, iarrayindex)
  read(svalue,*) ivalue

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_getvalue_int_direct (rparlist, ssectionName, &
                                         sparameter, ivalue, idefault,&
                                         iarrayindex)
!<description>
  
  ! Returns the value of a parameter in the section ssection. If the
  ! value does not exist, idefault is returned.  If idefault is not
  ! given, an error will be thrown.
  !
  ! If the value is an array of integers, the optional parameter
  ! iarrayindex>=0 allows to specify the number of the integer to be
  ! returned; iarrayindex=0 returns the value directly behind the '='
  ! sign in the line of the parameter, iarrayindex>0 returns the
  ! array-entry in the lines below the parameter.
  !
  ! When omitting iarrayindex, the value directly behind the '=' sign
  ! is returned.
  
!</description>
  
!<input>
    
  ! The parameter list.
  type(t_parlist), intent(IN) :: rparlist
  
  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! OPTIONAL: A default value
  integer, intent(IN), optional :: idefault
  
  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(IN), optional :: iarrayindex

!</input>
  
!<output>

  ! The value of the parameter
  integer, intent(OUT) :: ivalue
  
!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_MLDATA) :: sdefault,svalue
  
  ! Call the string routine, perform a conversion afterwards.
  if (present(idefault)) then
    write (sdefault,*) idefault
    call parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, sdefault, &
                                        iarrayindex)
  else
    call parlst_getvalue_string_direct (rparlist, ssectionName, &
                                        sparameter, svalue, &
                                        isubstring=iarrayindex)
  end if

  read(svalue,*) ivalue

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine parlst_addvalue_indir (rsection, sparameter, svalue, nsubstrings)
  
!<description>
  ! Adds a parameter to a section rsection.
  ! If the parameter exists, it's overwritten.
!</description>
  
!<inputoutput> 
    
  ! The section where to arr the parameter
  type(t_parlstSection), intent(INOUT) :: rsection
  
!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! The value of the parameter
  character(LEN=*), intent(IN) :: svalue
  
  ! OPTIONAL: Number of substrings. This allows a parameter to have
  ! multiple substrings, which can be accessed via the 'isubstring'
  ! parameter in the GET-routines.
  integer, intent(IN), optional :: nsubstrings
  
!</input>

!</subroutine>

    ! local variables
    character(LEN=PARLST_MLNAME) :: paramname
    integer :: i
    
    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call parlst_toupper (paramname)

    ! Get the parameter index into 'exists', finish.
    call parlst_fetchparameter(rsection, paramname, i)
    
    if (i .eq. 0) then
    
      ! Does not exist. Append.
      !
      ! Enough space free? Otherwise reallocate the parameter list
      if (rsection%iparamCount .eq. size(rsection%p_Sparameters)) then
        call parlst_reallocsection (rsection, size(rsection%p_Sparameters)+PARLST_NPARSPERBLOCK)
      end if

      ! Add the parameter - without any adjustment of the 'value' string
      rsection%iparamCount = rsection%iparamCount + 1  
      
      ! Set i to the index of the parameter
      i = rsection%iparamCount
      
    else
    
      ! Check if there are substrings. If yes, deallocate.
      ! Will be allocated later if necessary.
      if (associated(rsection%p_Rvalues(i)%p_Sentry)) then
        deallocate(rsection%p_Rvalues(i)%p_Sentry)
        rsection%p_Rvalues(i)%nsize = 0
      end if
      
    end if
      
    rsection%p_Sparameters(i) = paramname
    rsection%p_Rvalues(i)%sentry = svalue
    
    ! Add a list for the substrings if the parameter should have substrings.
    if (present(nsubstrings)) then
      if (nsubstrings .gt. 0) then
        allocate(rsection%p_Rvalues(i)%p_Sentry(nsubstrings))
        rsection%p_Rvalues(i)%nsize = nsubstrings
      else
        nullify(rsection%p_Rvalues(i)%p_Sentry)
        rsection%p_Rvalues(i)%nsize = 0
      end if
    else
      nullify(rsection%p_Rvalues(i)%p_Sentry)
      rsection%p_Rvalues(i)%nsize = 0
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_addvalue_direct (rparlist, ssectionName, sparameter, svalue,&
                                     nsubstrings)
!<description>
  
  ! Adds a parameter to a section with name ssectionName in the parameter list
  ! rparlist. If ssectionName='', the parameter is added to the unnamed
  ! section.
  
!</description>
  
!<inputoutput> 
    
  ! The parameter list.
  type(t_parlist), intent(INOUT) :: rparlist
  
!</inputoutput>

!<input>

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! The value of the parameter
  character(LEN=*), intent(IN) :: svalue
  
  ! OPTIONAL: Number of substrings. This allows a parameter to have
  ! multiple substrings, which can be accessed via the 'isubstring'
  ! parameter in the GET-routines.
  integer, intent(IN), optional :: nsubstrings

!</input>

!</subroutine>

  ! local variables
  type(t_parlstSection), pointer :: p_rsection
  
  ! Cancel if the list is not initialised.
  if (rparlist%isectionCount .eq. 0) then
    call output_line ('Parameter list not initialised!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_getvalue_string_fetch')
    call sys_halt()
  end if
  
  ! Get the section
  call parlst_querysection(rparlist, ssectionName, p_rsection) 
  if (.not. associated(p_rsection)) then
    call output_line ('Section not found: '//trim(ssectionName), &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_addvalue_direct')
    return
  end if

  ! Add the parameter 
  
  call parlst_addvalue_indir (p_rsection, sparameter, svalue, nsubstrings)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  subroutine parlst_setvalue_fetch (rsection, iparameter, svalue, iexists,&
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
  ! change. If omitted or = 0, the 'headline' directly behind the '='
  ! sign of the line 'name=value' is modified. Otherwise, the corresponding
  ! substring is changed.
  
!</description>
  
!<inputoutput> 
    
  ! The section where to arr the parameter
  type(t_parlstSection), intent(INOUT) :: rsection
  
!</inputoutput>

!<input>

  ! The parameter name.
  integer, intent(IN) :: iparameter

  ! The new value of the parameter
  character(LEN=*), intent(IN) :: svalue
  
  ! OPTIONAL: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  integer, intent(IN), optional :: isubstring

!</input>

!<output>

  ! Optional parameter. Is set to YES/NO, depending on whether
  ! the parameter exists.
  integer, intent(OUT), optional :: iexists

!</output>

!</subroutine>

  integer :: isub

  ! Check if iparameter is out of bounds. If yes, probably
  ! throw an error.
  
  if ((iparameter .lt. 0) .or. (iparameter .gt. rsection%iparamCount)) then
  
    if (.not. present(iexists)) then 
      call output_line ('Error. Parameter '//trim(sys_siL(iparameter,10))//&
          ' does not exist!', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_setvalue_fetch')
      call sys_halt()
    else
      iexists = NO
      return
    end if
  
  end if

  ! Depending on isubstring, change either the 'headline' or one
  ! of the substrings.
  isub = 0
  if (present(isubstring)) isub = isubstring

  if ((isub .le. 0) .or. &
      (isub .gt. rsection%p_Rvalues(iparameter)%nsize)) then
    rsection%p_Rvalues(iparameter)%sentry = svalue
  else
    rsection%p_Rvalues(iparameter)%p_Sentry(isub) = svalue
  end if

  if (present(iexists)) iexists = YES

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine parlst_setvalue_indir (rsection, sparameter, svalue, isubstring)
  
!<description>
  
  ! Modifies the value of a parameter in the section rsection.
  ! If the parameter does not exist, an error is thrown.
  !
  ! isubstring allows to specify the numer of a substring of the parameter to
  ! change. If omitted or = 0, the 'headline' directly behind the '='
  ! sign of the line 'name=value' is modified. Otherwise, the corresponding
  ! substring is changed.
  
!</description>
  
!<inputoutput> 
    
  ! The section where to arr the parameter
  type(t_parlstSection), intent(INOUT) :: rsection
  
!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! The new value of the parameter
  character(LEN=*), intent(IN) :: svalue
  
  ! OPTIONAL: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  integer, intent(IN), optional :: isubstring

!</input>

!</subroutine>

  ! local variables
  integer :: i,isub
  character(LEN=PARLST_MLNAME) :: paramname
  
  ! Create the upper-case parameter name
  paramname = adjustl(sparameter)
  call parlst_toupper (paramname)

  ! Get the parameter position
  i = parlst_queryvalue_indir (rsection, paramname)
  
  if (i .eq. 0) then
    call output_line ('Parameter '//trim(paramname)//&
        ' does not exist, cannot be modified!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_setvalue_indir')
    call sys_halt()
  else 
  
    ! Depending on isubstring, change either the 'headline' or one
    ! of the substrings.
    isub = 0
    if (present(isubstring)) isub = isubstring

    if ((isub .le. 0) .or. (isub .gt. rsection%p_Rvalues(i)%nsize)) then
      rsection%p_Rvalues(i)%sentry = svalue
    else
      rsection%p_Rvalues(i)%p_Sentry(isub) = svalue
    end if
  
  end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine parlst_setvalue_direct (rparlist, ssectionName, sparameter, svalue,&
                                     isubstring)
!<description>
  
  ! Modifies the value of a parameter in the section with name ssectionName
  ! in the parameter list rparlist.
  ! If the parameter does not exist, an error is thrown.
  !
  ! isubstring allows to specify the numer of a substring of the parameter to
  ! change. If omitted or = 0, the 'headline' directly behind the '='
  ! sign of the line 'name=value' is modified. Otherwise, the corresponding
  ! substring is changed.
  
!</description>
  
!<inputoutput> 
    
  ! The parameter list.
  type(t_parlist), intent(INOUT) :: rparlist
  
!</inputoutput>

!<input>

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(IN) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(IN) :: sparameter

  ! The new value of the parameter
  character(LEN=*), intent(IN) :: svalue
  
  ! OPTIONAL: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  integer, intent(IN), optional :: isubstring

!</input>

!</subroutine>

  ! local variables
  type(t_parlstSection), pointer :: p_rsection
  
  ! Cancel if the list is not initialised.
  if (rparlist%isectionCount .eq. 0) then
    call output_line ('Parameter list not initialised!', &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_setvalue_direct')
    call sys_halt()
  end if
  
  ! Get the section
  call parlst_querysection(rparlist, ssectionName, p_rsection) 
  if (.not. associated(p_rsection)) then
    call output_line ('Section not found: '//trim(ssectionName), &
        OU_CLASS_ERROR,OU_MODE_STD,'parlst_setvalue_direct')
    return
  end if

  ! Set the parameter 
  
  call parlst_setvalue_indir (p_rsection, sparameter, svalue, isubstring)

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Read a line from a text file.
  
  subroutine parlst_readlinefromfile (iunit, sdata, ilinelen, ios)
  
  ! The unit where to read from; must be connected to a file.
  integer, intent(IN) :: iunit
  
  ! The string where to write data to
  character(LEN=*), intent(OUT) :: sdata
  
  ! Length of the output
  integer, intent(OUT) :: ilinelen
  
  ! Status of the reading process. Set to a value <> 0 if the end
  ! of the file is reached.
  integer, intent(OUT) :: ios
  
  ! local variables
  integer :: eol
  character :: c
  
  sdata = ''
  ilinelen = 0
  
  ! Read the data - as long as the line/file does not end.
  eol = NO
  ios = 0
  do while ((ios .eq. 0) .and. (eol .eq. NO))
    
    ! Read a character.
    ! Unfortunately, Fortran forces me to use this dirty GOTO
    ! to decide processor-independently whether the line or
    ! the record ends.
    read (unit=iunit,fmt='(A1)',iostat=ios,advance='NO', end=10, eor=20) c
    goto 30
    
10  continue
    ! End of file. 
    ios = -1
    goto 30
    
20  continue
    ! End of record = END OF LINE.
    eol = YES

    ! Set error flag back to 0.
    ios = 0
    
30  continue    
    ! Don't do anything in case of an error
    if (ios .eq. 0) then
    
      ilinelen = ilinelen + 1
      sdata (ilinelen:ilinelen) = c
    
    end if
  
  end do
  
  end subroutine

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
  
  subroutine parlst_parseline (sdata, ityp, isubstring, ilinenum, &
                               ssecname, sparamname, svalue, sfilename)
  
  ! The line to be parsed
  character(LEN=*), intent(IN) :: sdata
  
  ! The typ of the line
  integer, intent(OUT) :: ityp

  ! input: =0: parse line as parameter. isubstring is changed to a value > 0
  !            is the parameter has multiple values attached.
  !        >0: parse line as substring of a multi-valued parameter, not 
  !            containing a leading 'name='.
  ! output: If the 'headline' of a multi-valued parameter is read, isubstring is
  !         changed to the number of substrings (the k in 'name(k)=...').
  !         Otherwise unchanged.
  integer, intent(INOUT) :: isubstring

  ! Line number
  integer, intent(IN) :: ilinenum
  
  ! Section name, if it's a section
  character(LEN=*), intent(INOUT) :: ssecname
  
  ! Parameter name, if it's a parameter
  character(LEN=*), intent(INOUT) :: sparamname
  
  ! Parameter value, if it's a parameter
  character(LEN=*), intent(INOUT) :: svalue
  
  ! OPTIONAL: Filename of the file to be parsed.
  ! Will be printed out in error messages if present.
  character(LEN=*), intent(IN), optional :: sfilename
  
  ! local variables
  integer :: i,j1,j2,ltr
  character(LEN=PARLST_LENLINEBUF) :: sbuf,slen
  
    ityp = 0
    
    ! Do we have data in sdata?
    if (sdata .eq. '') return
    
    ! Copy the input string - left adjusted - and get the string length
    sbuf = adjustl(sdata)
    
    ! Should we parse the line as first line of a parameter or as substring
    ! of a multi-valued parameter?
    if (isubstring .eq. 0) then
    
      ! Standard parameter or section header.
      !    
      ! Do we start with '[' and end with ']'?
      if (sbuf(1:1) .eq. "[") then
      
        ! Find the final ']'.
        do ltr = 1,len(sbuf)
          if (sbuf(ltr:ltr) .eq. "]") exit
        end do
        
        if (sbuf(ltr:ltr) .ne. ']') then
          if (present(sfilename)) then
            call output_line ('File: '//trim(sfilename),&
                OU_CLASS_ERROR,OU_MODE_STD,'parlst_parseline')
          end if
          call output_line ('Wrong syntax of section name. Line '//&
              trim(sys_siL(ilinenum,10))//':', &
              OU_CLASS_ERROR,OU_MODE_STD,'parlst_parseline')
          call output_line (sbuf, &
              OU_CLASS_ERROR,OU_MODE_STD,'parlst_parseline')
          call sys_halt()
        end if
        
        ! Get the section name
        ssecname = sbuf(2:ltr-1)
        ityp = 1
        return
        
      else if (sbuf(1:1) .eq. PARLST_COMMENT) then
      
        ! Comment sign
        return
        
      else
      
        ! Must be a parameter. Get the length of the string without comment
        ! at the end.
        call linelength(sbuf, ltr)
        
        ! ltr=0 means: empty line. Ignore that.
        if (ltr .eq. 0) return
        
        ! Is there a '(..)' that is indicating a multi-valued parameter?
        j1 = index(sbuf(1:ltr),'(')
        j2 = index(sbuf(1:ltr),')')

        ! Is there a '=' sign?
        i = index(sbuf(1:ltr),'=')

        if (i .eq. 0) then
          if (present(sfilename)) then
            call output_line ('File: '//trim(sfilename),&
                OU_CLASS_ERROR,OU_MODE_STD,'parlst_parseline')
          end if
          call output_line ('Invalid parameter syntax. Line '&
              //trim(sys_siL(ilinenum,10))//':', &
              OU_CLASS_ERROR,OU_MODE_STD,'parlst_parseline')
          call output_line (trim(sbuf), &
              OU_CLASS_ERROR,OU_MODE_STD,'parlst_parseline')
          call sys_halt()
        end if
      
        if ((j1 .eq. 0) .or. (j2 .le. j1)) then
        
          ityp = 2
          
          ! Get the name of the parameter
          sparamname = adjustl(sbuf(1:i-1))
          
          ! Get the parameter value
          svalue = adjustl(sbuf(i+1:ltr))
          
        else
        
          ! Probably multi-valued parameter with substrings in the 
          ! following lines.

          ! Get the name of the parameter
          sparamname = adjustl(sbuf(1:j1-1))

          ! Get the parameter value
          svalue = adjustl(sbuf(i+1:ltr))

          ! Get the length of the parameter list.
          slen = sbuf (j1+1:min(j2-1,len(slen)))

          isubstring = 0
          read(slen,*) isubstring
          
          if (isubstring .le. 0) then
            ! Oh, only one line. User want's to cheat :-)
            isubstring = 0
            
            ityp = 2
          else
            ! Real multi-valued parameter.
            ityp = 3
          end if
        
        end if
      
      end if
      
    else
      
      ! Substring of a multi-valued parameter.
      if (sbuf(1:1) .eq. PARLST_COMMENT) then
      
        ! Comment sign
        return
        
      else
       
        ! Must be a parameter. Get the length of the string without comment
        ! at the end.
        call linelength(sbuf, ltr)
        
        ! ltr=0 means: empty line. Ignore that.
        if (ltr .eq. 0) return
        
        ityp = 4
        
        ! Get the parameter value. Don't get a parameter name; there is none.
        svalue = adjustl(sbuf(1:ltr))
       
      end if
    
    end if
  
  contains
    
    ! Sub-subroutine: find the length of the line, removing comments
    ! at the end.
    
    subroutine linelength (sdata, l)
    
    ! The string to parse. Must not be ''!
    character(LEN=*), intent(IN) :: sdata
    
    ! The index of the last character without any comment at the end.
    integer, intent(OUT) :: l
    
    ! local variables
    logical :: bflag   ! Set to true if we are in apostroph mode
    integer :: lsdata
    
    bflag = .false.
    
    ! Go through all characters
    l = 0
    lsdata = len(sdata)
    do while (l .lt. lsdata)
      
      ! next character
      l = l+1

      ! A comment character while we are not in apostroph mode? Stop.
      if ((.not. bflag) .and. (sdata(l:l) .eq. PARLST_COMMENT)) then
        l = l-1
        exit
      end if
      
      ! An apostroph? 
      if (sdata(l:l) .eq. "'") then
      
        ! Switch the apostroph mode.
        ! Btw.: Two subsequent apostrophes will switch the mode off and on again.
        bflag = .not. bflag
      
      end if

    end do
    
    end subroutine
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine parlst_readfromfile (rparlist, sfilename, sdirectory)
  
!<description>
  
  ! This routine parses a text file for data of the INI-file form.
  ! sfilename must be the name of a file on the hard disc.
  ! The file may have references to subfiles in the unnamed section:
  ! If there is a parameter "simportdatafiles(n)" present in the
  ! unnamed section, the routine expects a list of filenames in this
  ! parameter. All the files listed there are read in as well.
  ! Parameters at the end of the master file will overwrite
  ! the parameters from the files in simportdatafiles.
  !
  ! The parameters read from the file(s) are added to the parameter list
  ! rparlist, which has to be initialised with parlst_init before
  ! calling the routine.
  !
  ! Remark: When adding parameters/sections to rparlist, the routine
  !   checks whether the parameters/sections already exist.
  !   Adding a parameter/section which exists does not result in an error -
  !   the first instance of the parameter/section will just be overwritten.
  
!</description>

!<inputoutput> 
    
  ! The parameter list which is filled with data from the file
  type(t_parlist), intent(INOUT) :: rparlist
  
!</inputoutput>
  
!<input>
  
  ! The filename of the file to read.
  character(LEN=*), intent(IN) :: sfilename
  
  ! OPTIONAL: Directory containing other data files in case
  ! the main file sfilename contains references to subfiles.
  character(LEN=*), intent(IN), optional :: sdirectory
  
!</input>

    ! local variables
    integer :: iidxsubfiles,icurrentsubfile
    character(LEN=PARLST_MLDATA), dimension(:), pointer :: p_Ssubfiles,p_SsubfilesTemp
    character(LEN=PARLST_MLDATA) :: sstring,smainfile
    integer :: nsubfiles,nnewsubfiles,j
    logical :: bexists

    ! Filename/path of the master dat file.
    ! Search at first in the specified path.
    bexists = .false.
    if (present(sdirectory)) then
      smainfile = trim(sdirectory)//"/"//sfilename
      inquire(file=smainfile, exist=bexists)
    end if
    
    if (.not. bexists) then
      smainfile = sfilename
      inquire(file=smainfile, exist=bexists)
    end if

    ! Create a list of files to be read.
    allocate(p_Ssubfiles(1))
    p_Ssubfiles(1) = smainfile
    
    icurrentsubfile = 0
    nsubfiles = 1

    ! Now read all files in the list. Append new files to the list if necessary.
    do while (icurrentsubfile .lt. nsubfiles)
      
      ! Read the unnamed section from the next file.
      icurrentsubfile = icurrentsubfile + 1

      ! Get the filename including the path. 
      bexists = .false.
      if (present(sdirectory)) then
        ! First search in the specified directory.
        sstring = trim(sdirectory)//"/"//trim(p_Ssubfiles(icurrentsubfile))
        inquire(file=sstring, exist=bexists)
      end if
      
      if (.not. bexists) then
        ! If no directory is specified or if no file was found there,
        ! directly search for the file without specifying the
        ! directory. So we will also find files where an absolute
        ! path is specified.
        sstring = trim(p_Ssubfiles(icurrentsubfile))
        inquire(file=sstring, exist=bexists)
      end if
      
      if (bexists) then
        call parlst_readfromsinglefile (rparlist, sstring, .false.)
        
        ! Replace the filename with the string including the path.
        ! Then the actual read process at the end of the routine can be
        ! handled easier.
        p_Ssubfiles(icurrentsubfile) = trim(sstring)
        
      else
        call output_line ('Specified data-subfile does not exist: '//sstring, &
            OU_CLASS_WARNING,OU_MODE_STD,'parlst_readfromfile')
      end if

      ! Check if there's a parameter "simportdatafiles" available in the
      ! parameter list. It must be present in the unnamed section.
      call parlst_fetchparameter(rparlist%p_Rsections(1), "SIMPORTDATAFILES", iidxsubfiles)
      
      if (iidxsubfiles .ne. 0) then
        ! Append the new files to the file list.
        !
        ! Get the number of new files. The parameter definitely exists as
        ! it was created when reading the 'master' INI file.
        nnewsubfiles = parlst_querysubstrings (rparlist%p_Rsections(1), &
            "SIMPORTDATAFILES")
        
        ! if nnewsubfiles=0, there is (hopefully) only one string here.
        if (nnewsubfiles .eq. 0) then
          call parlst_getvalue_string(rparlist%p_Rsections(1), iidxsubfiles, sstring)
          if (trim(sstring) .ne. "") then
            ! Append the data.
            allocate(p_SsubfilesTemp(nsubfiles+1))
            p_SsubfilesTemp(1:nsubfiles) = p_Ssubfiles(:)
            deallocate(p_Ssubfiles)
            p_Ssubfiles => p_SsubfilesTemp
            nullify(p_SsubfilesTemp)
            
            read (sstring,*) p_Ssubfiles(nsubfiles+1)
            nsubfiles = nsubfiles + 1
          end if
        else
          ! Get all the filenames.
          allocate(p_SsubfilesTemp(nsubfiles+nnewsubfiles))
          p_SsubfilesTemp(1:nsubfiles) = p_Ssubfiles(:)
          deallocate(p_Ssubfiles)
          p_Ssubfiles => p_SsubfilesTemp
          nullify(p_SsubfilesTemp)
          
          do j=1,nnewsubfiles
            call parlst_getvalue_string(rparlist%p_Rsections(1), iidxsubfiles, &
                sstring,isubstring=j)
            read (sstring,*) p_Ssubfiles(nsubfiles+j)
          end do
          
          nsubfiles = nsubfiles + nnewsubfiles
        end if
        
        ! Remove all substrings from the simportdatafiles parameter, so the parameter
        ! is filled with new data upon the next read statement.
        call parlst_addvalue(rparlist%p_Rsections(1), "SIMPORTDATAFILES", "")
          
      end if
      
    end do
    
    ! Ok, at that point we know which files to read -- so read them, one after
    ! the other. The 'master' file must be read at last!
    do icurrentsubfile = 2,nsubfiles
      sstring = trim(p_Ssubfiles(icurrentsubfile))
      inquire(file=sstring, exist=bexists)
      
      if (bexists) then
        ! Read the sub-files...
        call parlst_readfromsinglefile (rparlist, sstring, .true.)
      end if
    end do
    
    ! ... and the master
    icurrentsubfile = 1
    sstring = trim(p_Ssubfiles(icurrentsubfile))
    inquire(file=sstring, exist=bexists)
    
    if (bexists) then
      call parlst_readfromsinglefile (rparlist, sstring, .true.)
    end if

    if (nsubfiles .gt. 1) then
      ! There have a couple of subfiles been read from disc.
      ! All subfiles can be found in p_Ssubfiles.
      
      ! Incorporate the complete list to the parameter "SIMPORTDATAFILES".
      if (associated(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_Sentry)) then
        deallocate(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_Sentry)
      end if
      
      allocate(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_Sentry(nsubfiles-1))
      do icurrentsubfile = 1,nsubfiles-1
        rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_Sentry(icurrentsubfile) = &
          '"'//trim(p_Ssubfiles(1+icurrentsubfile))//'"'
      end do
      rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%nsize = nsubfiles-1
    end if

    ! Release memory, finish
    deallocate(p_Ssubfiles)
    
  end subroutine
  
  ! ***************************************************************************

  subroutine parlst_readfromsinglefile (rparlist, sfilename, bimportSections)
    
!<description>
  
  ! This routine parses a text file for data of the INI-file form.
  ! sfilename must be the name of a file on the hard disc.
  ! The parameters read from the file are added to the parameter list
  ! rparlist, which has to be initialised with parlst_init before
  ! calling the routine.
  ! Remark: When adding parameters/sections to rparlist, the routine
  !   checks whether the parameters/sections already exist.
  !   Adding a parameter/section which exists does not result in an error -
  !   the first instance of the parameter/section will just be overwritten.
  
!</description>

!<inputoutput> 
    
  ! The parameter list which is filled with data from the file
  type(t_parlist), intent(INOUT) :: rparlist
  
!</inputoutput>
  
!<input>
  
  ! The filename of the file to read.
  character(LEN=*), intent(IN) :: sfilename
  
  ! TRUE: Import all sections in the DAT file.
  ! FALSE: Import only the main (unnamed) section and ignore all other
  ! sections.#
  logical, intent(in) :: bimportSections
  
!</input>
    
    ! local variables
    integer :: iunit,ios,isbuflen,ityp,ilinenum,isubstring,nsubstrings,iparpos
    type(t_parlstSection), pointer :: p_currentsection
    character(LEN=PARLST_LENLINEBUF) :: sdata
    character(LEN=PARLST_MLSECTION) :: ssectionname
    character(LEN=PARLST_MLNAME) :: sparname
    character(LEN=PARLST_MLDATA) :: svalue
    
    ! Try to open the file
    call io_openFileForReading(sfilename, iunit)
    
    ! Oops...
    if (iunit .eq. -1) then
      call output_line ('Error opening .INI file: '//trim(sfilename), &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_readfromsinglefile')
      call sys_halt()
    end if
    
    ! Start adding parameters to the unnamed section
    p_currentsection => rparlist%p_Rsections(1)
    
    ! Read all lines from the file
    ios = 0
    ilinenum = 0
    isubstring = 0
    nsubstrings = 0
    do while (ios .eq. 0) 
      
      ! Read a line from the file into sbuf
      call parlst_readlinefromfile (iunit, sdata, isbuflen, ios)
      ilinenum = ilinenum + 1
      
      if (isbuflen .ne. 0) then
      
        ! Parse the line
        call parlst_parseline (sdata, ityp, nsubstrings, ilinenum, ssectionname, &
                              sparname, svalue,sfilename)  
        
        select case (ityp)
        case (1)
          ! Stop parsing the file here if bimportSections tells us to do so.
          ! When the first section starts, the unnamed section is finished.
          if (.not. bimportSections) exit
        
          ! Check if the section exists; if not, create a new one.
          call parlst_querysection(rparlist, ssectionname, p_currentsection) 
          
          if (.not. associated(p_currentsection)) then
            ! A new section name. Add a section, set the current section
            ! to the new one.
            call parlst_addsection (rparlist, ssectionname)
            p_currentsection => rparlist%p_Rsections(rparlist%isectionCount)
          end if
          
        case (2)
          ! A new parameter. Add it to the current section.
          call parlst_addvalue (p_currentsection, sparname, svalue)
          
        case (3)
          ! 'Headline' of a multi-valued parameter. Add the parameter with
          ! isubstring subvalues
          call parlst_addvalue (p_currentsection, sparname, svalue, nsubstrings)
          
          ! Fetch the parameter for later adding of subvalues.
          iparpos = parlst_queryvalue(p_currentsection, sparname)
          
          ! isubstring counts the current readed substring.
          ! Set it to 0, it will be increased up to nsubstrings in 'case 4'.
          isubstring = 0
          
        case (4)
          ! Increase number of current substring
          isubstring = isubstring + 1
          
          ! Sub-parameter of a multi-valued parameter. Add the value to
          ! the last parameter that was added in case 3.
          call parlst_setvalue_fetch (p_currentsection, iparpos, svalue, &
                                      isubstring=isubstring)
                                      
          ! Decrement the substring counter. If we reach 0, parlst_parseline
          ! continues to parse standard parameters.
          nsubstrings = nsubstrings - 1
          
        ! Other cases: comment.
        end select
      
      end if
    
    end do
    
    ! Close the file.
    close (iunit)

    ! Expand all subvariables to the actual values.
    call parlst_expandsubvars(rparlist)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine parlst_getStringRepresentation (rparlist, p_sconfiguration)
  
!<description>
  ! Creates a string representation of the given parameter list rparlist.
  ! p_sconfiguration will be created as "array[1..*] of char" on
  ! the heap containing this representation. The memory must be manually 
  ! released by the caller using DEALLOCATE when finished using the string
  ! representation.
!</description>
  
!<input> 
  ! The parameter list which is filled with data from the file
  type(t_parlist), intent(IN) :: rparlist
!</input>

!<output>
  ! A pointer to a character array containing all lines of the parameter list.
  ! Points to NULL() if there is no data in the parameter list.
  ! Each line is terminated by NEWLINE.
  ! If there is data, a new pointer is allocated for this on the heap. 
  ! The user must manually release the memory when finished using it.
  character, dimension(:), pointer :: p_sconfiguration
!</output>

!</subroutine>

  integer :: ilength,isection,ivalue,ientry,icount
  character, dimension(:), pointer :: p_sbuf
  
    if (rparlist%isectionCount .eq. 0) then
      call output_line ('Parameter list not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_getStringRepresentation')
      call sys_halt()
    end if
  
    nullify(p_sbuf)
    
    ! Number of characters in the buffer
    ilength = 0
    
    ! Loop through all sections
    do isection = 1,rparlist%isectionCount
    
      ! Append the section name. May be empty for the unnamed section,
      ! which is always the first one.
      if (isection .gt. 1) then
        ! Empty line before
        if (ilength .gt. 0) call appendString(p_sbuf,ilength,'')
        call appendString(p_sbuf,ilength,&
          '['//trim(rparlist%p_Rsections(isection)%ssectionName)//']')
      end if
        
      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount
        
        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        if (icount .eq. 0) then
          ! Write "name=value"
          call appendString(p_sbuf,ilength,&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"="// &
            trim(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%sentry))
        else
          ! Write "name(icount)="
          call appendString(p_sbuf,ilength,&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//trim(sys_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          do ientry = 1,icount
            call appendString(p_sbuf,ilength,&
              trim(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
                   p_Sentry(ientry)))
          end do
        end if
      
      end do ! ivalue
    
    end do ! isection
    
    ! Allocate a new character array with the correct size, copy p_sbuf to
    ! that ald release the old p_sbuf.
    ! Return NULL() if there is no data.
    nullify(p_sconfiguration)
    if (ilength .gt. 0) then
      allocate(p_sconfiguration(ilength))
      p_sconfiguration = p_sbuf(1:ilength)
    end if
    
    ! Release our temp buffer
    if (associated(p_sbuf)) deallocate(p_sbuf)

  contains
  
    ! Makes sure, the character buffer points to a character memory block of
    ! size nsize. If not, the block is reallocated to have that size.
    subroutine assumeBufSize(p_sconfig,nsize)
    
    character, dimension(:), pointer :: p_sconfig
    integer, intent(IN) :: nsize
    
    character, dimension(:), pointer :: p_sconfignew
    
      if (.not. associated(p_sconfig)) then
        allocate(p_sconfig(nsize))
      else if (size(p_sconfig) .lt. nsize) then
        allocate(p_sconfignew(nsize))
        p_sconfignew(1:size(p_sconfig)) = p_sconfig
        deallocate(p_sconfig)
        p_sconfig => p_sconfignew
      end if
    
    end subroutine
    
    ! Appends sstring to the buffer p_sconfig, followed by a NEWLINE
    ! character. Reallocates memory if necessary.
    subroutine appendString(p_sconfig,iconfigLength,sstring)
    
    ! Pointer to character data
    character, dimension(:), pointer :: p_sconfig
    
    ! In: Current length of data stream in p_sconfig.
    ! Out: New length of data stream in p_sconfig
    integer, intent(INOUT) :: iconfigLength
    
    ! The string to be added.
    character(LEN=*), intent(IN) :: sstring
    
      integer :: nblocks,nblocksneeded,i
    
      ! How many memory blocks do we need for the current configuration?
      ! We work block-wise to prevent too often reallocation.
      if (.not. associated(p_sconfig)) then
        nblocks = 0
      else
        nblocks = size(p_sconfig) / SYS_STRLEN
      end if
      nblocksneeded = 1 + (iconfigLength+len(sstring)+1) / SYS_STRLEN
      if (nblocksneeded .gt. nblocks) then
        call assumeBufSize(p_sconfig,nblocksneeded*SYS_STRLEN)
      end if
      
      ! Append the data
      do i=1,len(sstring)
        iconfigLength = iconfigLength+1
        p_sconfig(iconfigLength) = sstring(i:i)
      end do
      
      ! Append NEWLINE as line-end character
      iconfigLength = iconfigLength+1
      p_sconfig(iconfigLength) = NEWLINE
    
    end subroutine
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine parlst_info (rparlist)
  
!<description>
  ! Prints the parameter list rparlist to the terminal.
!</description>
  
!<input> 
  ! The parameter list which is to be printed to the terminal.
  type(t_parlist), intent(IN) :: rparlist
!</input>

!</subroutine>


  integer :: isection,ivalue,ientry,icount
  
    if (rparlist%isectionCount .eq. 0) then
      call output_line ('Parameter list not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_info')
      call sys_halt()
    end if
  
    ! Loop through all sections
    do isection = 1,rparlist%isectionCount
    
      ! Append the section name. May be empty for the unnamed section,
      ! which is always the first one.
      if (isection .gt. 1) then
        ! Empty line before
        call output_lbrk()
        call output_line('['//trim(rparlist%p_Rsections(isection)%ssectionName)//']')
      end if
        
      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount
        
        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        if (icount .eq. 0) then
          ! Write "name=value"
          call output_line(&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"="// &
            trim(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%sentry))
        else
          ! Write "name(icount)="
          call output_line(&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//trim(sys_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          do ientry = 1,icount
            call output_line(&
              trim(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
                   p_Sentry(ientry)))
          end do
        end if
      
      end do ! ivalue
    
    end do ! isection
    
  end subroutine

  ! ***************************************************************************

!<function>
  function parlst_expandEnvVariables(sbuffer) result(sresult)

  !<description>
    ! This subroutine recursively expands all environment variables in the given
    ! string sbuffer.
  !</description>

  !<input>
    ! string
    character(len=*), intent(in) :: sbuffer
  !</input>

  !<result>
    character(len=SYS_STRLEN) :: sresult
  !</result>

  !<errors>
    ! ERR_CONV_ENV_VAR_NOT_FOUND (critical)
  !</errors>
!</function>

    ! flag
    logical :: bfoundInEnv

    ! start and end position of variable
    integer(I32) :: istartPos, istopPosRelative

    ! variable to expand environment variable on-the-fly to if found
    character(len=SYS_STRLEN) :: sauxEnv


    ! Initialise return value
    sresult = trim(sbuffer)

    ! check for a $ character
    istartPos = index(sresult, "$")
    do while (istartPos .gt. 0)
      ! Detect end of variable: a variable ends at the first character that
      ! is neither in '[A-Z]', '[0-9]' nor '_'.
      istopPosRelative = verify(sresult(istartPos+1:), &
                          "abcdefghijklmnopqrstuvwxyz" // &
                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ" // &
                          "0123456789_")

      bfoundInEnv = .FALSE.
      ! Retrieve value of environment variable
      ! (Do not forget to cut the dollar sign.)
      bfoundInEnv = &
           parlst_getenv_string(trim(&
               sresult(istartPos + 1 : istartPos + istopPosRelative - 1)), sauxEnv)
      if (bfoundInEnv) then
        ! Replace environment variable by its content
        sresult = sresult(1:istartPos-1) // &
                  trim(sauxEnv) // &
                  trim(sresult(istartPos + istopPosRelative:))
      else
        call error_print(ERR_PARLST_ENV_VAR_NOT_FOUND, "parlst_expandEnvVariables", &
                         ERR_CRITICAL, &
                         sarg1=trim(sresult(istartPos + 1 : istartPos + istopPosRelative - 1)))
      endif

      ! check for next $ character
      istartPos = index(sresult, "$")
    enddo

  end function parlst_expandEnvVariables

  ! ***************************************************************************

!<function>
  logical function parlst_getenv_string(svar, sresult)

  !<description>
    ! This functions returns the string value of a given enviroment variable. The routine
    ! returns .TRUE., if the variable exists, otherwise .FALSE. .
  !</description>

  !<input>
    ! name of the enviroment variable
    character(len=*), intent(in) :: svar
  !</input>

  !<output>
    ! value of the enviroment variable
    character(len=SYS_STRLEN), intent(out) :: sresult
  !</output>

  !<result>
    ! exit status
  !</result>

  !<errors>
    ! none
  !</errors>
!</function>

    character(len=SYS_STRLEN) :: svalueInEnv

    integer :: nstatus

#if defined USE_COMPILER_NEC || defined USE_COMPILER_PGI_6_1 || defined USE_COMPILER_PGI_6_2 || defined USE_COMPILER_PGI_7_0
    call getenv(trim(svar), svalueInEnv)
    if (trim(svalueInEnv) .eq. "") then
      nstatus = 1
    else
      nstatus = 0
    endif
#else
    call get_environment_variable(trim(svar), svalueInEnv, status=nstatus)
#endif

    select case (nstatus)
    case (0)
      ! Copy string only up to first whitespace character
!      read(svalueInEnv, '(A)') sresult

      ! Copy complete string
      sresult = svalueInEnv
      parlst_getenv_string = .TRUE.

    case (1)
      ! Environment variable does not exist
      sresult = ""
      parlst_getenv_string = .FALSE.

    case default
      !  2: Processor does not support environment variables
      ! >2: Some error occurred
      ! -1: variable svalueInEnv too short to absorb environment variables` content
      sresult = ""
      parlst_getenv_string = .FALSE.

    end select

  end function parlst_getenv_string
  
  ! ***************************************************************************

!<subroutine>

  subroutine parlst_expandsubvars(rparlist)
  
!<description>
  ! This subroutine expands variables when referring to subvariables:
  ! The value of a parameter may refer to another variable in the parameter
  ! list. This can be specified by tokens of the form "%{NAME}" or "{SECTION.NAME}",
  ! depending on whether it is part of the main section or not.
  ! Example: "NLMIN = %{NLMAX}"; in this case, NLMIN is always set to
  ! the same value as NLMAX.
  !
  ! The routine parses all variables in the DAT file to resolve such
  ! references. Note that recursive definitions are not allowed!
!</description>
  
!<inputoutput> 
  ! The parameter list which is filled with data from the file
  type(t_parlist), intent(INOUT) :: rparlist
!</inputoutput>

!</subroutine>

    integer :: isection,ivalue,ientry,icount
    
    if (rparlist%isectionCount .eq. 0) then
      call output_line ('Parameter list not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_expandsubvars')
      call sys_halt()
    end if
  
    ! Loop through all sections
    do isection = 1,rparlist%isectionCount
    
      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount
        
        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        if (icount .eq. 0) then
          ! Expand the value if is refers to subvalues.
          call expandSubvariable(rparlist,&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%sentry)
        else
          ! Loop through the subvalues.
          do ientry = 1,icount
            ! Expand the value if is refers to subvalues.
            call expandSubvariable(rparlist,&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_Sentry(ientry))
          end do
        end if
      
      end do ! ivalue
    
    end do ! isection
      
  contains
  
    subroutine findSubvariable(sstring,istart,iend,ssection,sname,ivalue)
    
    ! Searchs in the string sstring for the first occurance of a subvariable.
    ! A subvariable has the form "%{NAME}" or "{SECTION.NAME}"
    ! or "%{NAME:INDEX}" or "%{SECTION.NAME:INDEX}"
    ! depending on whether it is part of the main section or not.
    
    ! The string where a subvariable is searched.
    character(len=*), intent(in) :: sstring
    
    ! Returns the start of the variable in the string or 0 if no subvariable
    ! is found.
    integer, intent(out) :: istart

    ! Returns the end of the variable in the string or 0 if no subvariable
    ! is found.
    integer, intent(out) :: iend
    
    ! Returns the name of the section or "" if either no subvariable is found
    ! or the unnamed section is referred to.
    character(len=*), intent(out) :: ssection

    ! Returns the name of the subvariable or "" if no subvariable is found.
    character(len=*), intent(out) :: sname
    
    ! Returns the number/index INDEX of the subvalue or 0, if there is
    ! no index or if no subvariable is found.
    integer, intent(out) :: ivalue
    
    ! local variables
    integer :: i,j,istrlen,idotpos,icolonpos
    logical :: bstartfound
    
      ! Ok, this is a parser. We have the following rules:
      ! "%%" means one "%" and is not interpreted as the begin of a
      ! token.
      ! "%{NAME}" is a token referring to a variable in the unnamed section.
      ! "%{NAME:INDEX}" is a token referring to a subvariable
      ! of variable NAME in section SECTION.
      ! "%{SECTION.NAME}" is a token referring to a name in a named section
      ! which must exist.
      ! "%{SECTION.NAME:INDEX}" is a token referring to a subvariable
      ! of variable NAME in section SECTION.
      !
      istart = 0
      iend = 0
      ivalue = 0
      bstartfound = .false.
      
      ! Lets loop through the characters.
      istrlen = LEN(sstring)
      do i=1,istrlen
        if (sstring(i:i) .eq. "%") then
          ! Did we already found the "%"?
          if (bstartfound) then
            ! That's our escape sequence. Don't do anything, just
            ! return to 'normal' mode.
            bstartfound = .false.
          else
            ! That's probably the beginning of a token.
            bstartfound = .true.
          end if
          
          ! Go on.
          cycle
        end if
        
        ! The next things only execute if we are close to a token...
        if (bstartfound) then
          if (sstring(i:I) .eq. "{") then
            ! Yes, that's a token.
            istart = i-1
            
            ! Find the end of the token and probably the dot/colon
            idotpos = 0
            icolonpos = 0
            do j=istart+1,istrlen
              if (sstring(j:J) .eq. ".") then
                ! Here is the dot.
                idotpos = j
              end if

              if (sstring(j:J) .eq. ":") then
                ! Here is the dot.
                icolonpos = j
              end if
              
              if (sstring(j:j) .eq. "}") then
                ! Here, the token ends. 
                iend = j
                
                ! Extract name and probably the section
                if (idotpos .eq. 0) then
                  ssection = ""
                  if (icolonpos .eq. 0) then
                    sname = sstring(istart+2:iend-1)
                  else
                    sname = sstring(istart+2:icolonpos-1)
                    
                    ! Get the value
                    read(sstring(icolonpos+1:iend-1),*) ivalue
                  end if
                else
                  ssection = sstring(istart+2:idotpos-1)
                  if (icolonpos .eq. 0) then
                    sname = sstring(idotpos+1:iend-1)
                  else
                    sname = sstring(idotpos+1:icolonpos-1)

                    ! Get the value
                    read(sstring(icolonpos+1:iend-1),*) ivalue
                  end if
                end if
                
                ! That's it.
                return
                
              end if
            end do
          end if
        end if
      end do
      
      ! Nothing found.
      ssection = ""
      sname = ""
      
    end subroutine
    
    ! ---------------------------------------------------------------
    
    subroutine expandSubvariable(rparlist,sstring)
    
    ! Expands the subvariables in sstring to a fully qualified string.
    
    ! The parameter list containing the variables that can be used
    ! as subvariables.
    type(t_parlist), intent(in) :: rparlist
    
    ! The string to be expanded. Receives the expanded string upon return.
    character(len=*), intent(inout) :: sstring
    
      ! local variables
      integer :: istart,iend,ilen,ivalue
      character(len=PARLST_MLSECTION) :: ssection
      character(len=PARLST_MLNAME) :: sname
      character(len=LEN(sstring)) :: sbuffer
      character(len=PARLST_MLDATA) :: sdata
      
      ! Repeat until we found all subvariables
      istart = 1
      iend = 0
      do while (istart .ne. 0)
  
        ! Find the first subvariable
        call findSubvariable(sstring,istart,iend,ssection,sname,ivalue)
        
        if (istart .ne. 0) then
          ! Copy the string to the buffer
          sbuffer = sstring
          
          ! Now copy back and replace the variable by the stuff from the
          ! parameter list.
          if (ivalue .eq. 0) then
            call parlst_getvalue_string (rparlist, ssection, sname, sdata)
          else
            call parlst_getvalue_string (rparlist, ssection, sname, sdata, &
                isubstring=ivalue)
          end if
          sstring = sbuffer(1:istart-1)//TRIM(sdata)//sbuffer(iend+1:)
          
        end if
      
      end do
    
    end subroutine

  end subroutine
  
end module
