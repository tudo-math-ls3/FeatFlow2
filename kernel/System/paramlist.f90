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
!# <verb>
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
!# </verb>
!# the .INI file is build up by different sections.
!# Each section starts with the section name enclosed by
!# brackets ('[...]'), followed by a list of parameters
!# if the form "name=parameter". There is one unnamed parameter
!# block at the beginning of the file.
!#
!# There is at most one parameter per line. The parameter names
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
!# The subvariables feature \\
!# ======================== \\
!# The parameter list also allows to specify variables as subvariables
!# of other variables. Take a look at the following example:
!#
!# <verb>
!# -------------------------snip------------------------------
!# NLMIN=1
!# ACOMPLEXSTRING=! %{NLMIN} ! %{imainelement} ! %{SECTION1.ilist:3}
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
!# </verb>
!# When reading the file, the parameter "ACOMPLEXSTRING" is automatically
!# expanded to the value "! 1 ! 5 ! ghi" by using other variables from
!# the parameter file. The following syntax is allowed here to refer
!# to other parameters:
!#
!# <verb>
!#  %{NAME}               - A variable from the unnamed section
!#  %{NAME:idx}           - Value number idx of variable NAME
!#                          from the unnamed section
!#  %{SECTION.NAME}       - Variable NAME from section SECTION
!#  %{SECTION.NAME:idx}   - Value number idx of variable NAME
!#                          from section SECTION
!# </verb>
!#
!# The environment variable feature \\
!# ================================ \\
!# DAT files may refer to environment variables. Example:
!# <verb>
!# -------------------------snip------------------------------
!# NLMIN=1
!# NLMAX=$NLMAXENV
!# -------------------------snip------------------------------
!# </verb>
!# In this case, it is assumed that $NLMAXENV is an environment
!# variable of the system where the program runs. So when NLMAX
!# is requested, the value of the environment variable $NLMAXENV
!# is returned.
!#
!# The subfile feature \\
!# =================== \\
!# An INI file may contain references to subfiles. Subfiles must
!# be specified at the beginning of an INI file with the following
!# syntax:
!# <verb>
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
!# </verb>
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
!# 14.) parlst_findvalue
!#      -> Checks if a given value exists and returns its substring index
!#
!# Auxiliary routines
!#
!# 1.) parlst_readfromsinglefile
!#     -> Reads a single INI file without checking for subfiles.
!#
!# 2.) parlst_expandEnvVariables
!#     -> expands all environment variables in the parameter list
!#
!# 3.) parlst_expandSubvars
!#     -> expands all subvariables in the parameter list
!#
!# </purpose>
!##############################################################################

module paramlist

!$ use omp_lib
  use fsystem
  use io
  use genoutput

  implicit none

  private

!<constants>

  !<constantblock>

  ! Maximum length of a section name.
  integer, parameter, public :: PARLST_MLSECTION = 64

  ! Maximum length of parameter names: 32 characters
  integer, parameter, public :: PARLST_MLNAME = 32

  ! Default length of parameter data: 256 characters
  integer, parameter, public :: PARLST_MLDATA = 256

  ! Minimum number of free parameter 'slots' per parameter section.
  ! If there are too many parameters in a parameter section, the
  ! structure is dynamically extended in terms of PARLST_NPARSPERBLOCK
  ! entries.
  integer, parameter, public :: PARLST_NPARSPERBLOCK = 32

  ! Minimum number of parameter sections.
  ! If there are too many parameter sections in a parameter block, the
  ! structure is dynamically extended in terms of PARLST_NSECTIONS
  ! entries.
  integer, parameter, public :: PARLST_NSECTIONS = 8

  ! Maximum length of a line in a INI file. Lines longer than this
  ! are truncated.
  integer, parameter, public :: PARLST_LENLINEBUF = 1024

  ! Comment character
  character, parameter, public :: PARLST_COMMENT = "#"

  !</constantblock>

!</constants>

!<types>

  !<typeblock>

  ! This structure realises a value associated to a parameter name.
  ! A value consists of one string or an array of strings.

  type t_parlstValue

    private

    ! Number of strings. If set to 0, the value consists of one
    ! string, to be found in svalue. If > 0, there are nsize
    ! strings to be found in p_SentryList.
    integer :: nsize = 0

    ! Single string; contains the value in case nsize=0
    character, dimension(:), pointer :: p_sentry => null()

    ! Array of strings in case nsize>0
    character, dimension(:,:), pointer :: p_SentryList => null()

  end type

  public :: t_parlstValue

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

  public :: t_parlstSection

  !</typeblock>

  !<typeblock>

  ! This structure realises a parameter list. Parameters can be read into
  ! it from a file. Parameters can be obtained from the structure using
  ! the query/get routines.

  type t_parlist

    private

    ! Actual number of sections in the parameter list. There is at least
    ! one section - the unnamed section. If this value is =0, the parameter
    ! list is not initialised.
    integer :: isectionCount = 0

    ! A list of sections. The first section is always the unnamed section.
    type(t_parlstSection), dimension(:), pointer :: p_Rsections => null()

  end type

  public :: t_parlist

  !</typeblock>

!</types>

  private :: parlst_initsection, parlst_reallocsection, parlst_realloclist
  private :: parlst_reallocSubVariables
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
    module procedure parlst_getvalue_int8_fetch
    module procedure parlst_getvalue_int16_fetch
    module procedure parlst_getvalue_int32_fetch
    module procedure parlst_getvalue_int64_fetch
    module procedure parlst_getvalue_int8_indir
    module procedure parlst_getvalue_int16_indir
    module procedure parlst_getvalue_int32_indir
    module procedure parlst_getvalue_int64_indir
    module procedure parlst_getvalue_int8_direct
    module procedure parlst_getvalue_int16_direct
    module procedure parlst_getvalue_int32_direct
    module procedure parlst_getvalue_int64_direct
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

  interface parlst_findvalue
    module procedure parlst_findvalue_direct
    module procedure parlst_findvalue_indir
  end interface

  public :: parlst_init
  public :: parlst_readfromfile
  public :: parlst_clear
  public :: parlst_done
  public :: parlst_querysection
  public :: parlst_addsection
  public :: parlst_queryvalue
  public :: parlst_querysubstrings
  public :: parlst_getvalue_string
  public :: parlst_getvalue_int
  public :: parlst_getvalue_single
  public :: parlst_getvalue_double
  public :: parlst_addvalue
  public :: parlst_setvalue
  public :: parlst_getStringRepresentation
  public :: parlst_info
  public :: parlst_readfromsinglefile
  public :: parlst_expandEnvVariables
  public :: parlst_expandSubvars
  public :: parlst_dumptofile
  public :: parlst_findvalue

contains

  ! ***************************************************************************

  ! Internal subroutine: Initialise a newly created parameter section.

  subroutine parlst_initsection (rparlstSection,sname)

  type(t_parlstSection), intent(inout) :: rparlstSection
  character(LEN=*), intent(in) :: sname

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
  type(t_parlstSection), intent(inout) :: rparlstSection

  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the section should be able to handle.
  integer, intent(in) :: inewsize

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

  ! Internal subroutine: Reallocates a sub-parameter list
  ! This increases the size of the sub-parameters of a parameter.
  ! Old strings are copied.

  subroutine parlst_reallocSubVariables (rvalue, inewsize)

  ! THe parameter item where to reallocate sub-parameters.
  type(t_parlstValue), intent(inout) :: rvalue

  ! The new 'length' of the sub-parameters.
  integer, intent(in) :: inewsize

    ! local variables
    integer :: i,iold
    character, dimension(:,:), pointer :: p_Sdata

    if (ubound(rvalue%p_SentryList,2) .eq. 0) return ! nothing to do

    ! Old size
    iold = ubound(rvalue%p_SentryList,1)

    ! Allocate memory for the strings.
    allocate(p_Sdata(inewsize,ubound(rvalue%p_SentryList,2)))
    p_Sdata(:,:) = ' '

    ! Copy the substrings.
    do i=1,ubound(rvalue%p_SentryList,2)
      p_Sdata(1:min(inewsize,iold),i) = rvalue%p_SentryList(1:min(inewsize,iold),i)
    end do

    ! Replace the data.
    deallocate(rvalue%p_SentryList)
    rvalue%p_SentryList => p_Sdata

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Release a section.
  ! Removes all temporary memory that is allocated by a section.

  subroutine parlst_releasesection (rparlstSection)

  ! The section to release.
  type(t_parlstSection), intent(inout) :: rparlstSection

  ! local variables
  integer :: i

  ! Loop through all values in the current section if there is
  ! an array-value. Release them.
  do i=size(rparlstSection%p_Rvalues),1,-1
    if (associated(rparlstSection%p_Rvalues(i)%p_sentry)) then
      deallocate(rparlstSection%p_Rvalues(i)%p_sentry)
    end if
    if (rparlstSection%p_Rvalues(i)%nsize .gt. 0) then
      deallocate(rparlstSection%p_Rvalues(i)%p_SentryList)
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
  type(t_parlist), intent(inout) :: rparlist

  ! The new 'size' of the section, i.e. the new number of parameters,
  ! the section should be able to handle.
  integer, intent(in) :: inewsize

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

  ! Internal subroutine: Search in a section for a parameter
  ! and return the index - or 0 if the parameter does not exist.

  subroutine parlst_fetchparameter(rsection, sname, iparamnum)

  ! The section.
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name to look for. Must be uppercase.
  character(LEN=*), intent(in) :: sname

  ! The number of the parameter in the list or 0 if it does not exist.
  integer, intent(out) :: iparamnum

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
  type(t_parlist), intent(inout) :: rparlist

!</inputoutput>

!</subroutine>

  ! Set the section-count to 1.
  rparlist%isectionCount = 1

  ! Allocate a first set of sections
  allocate(rparlist%p_Rsections(PARLST_NSECTIONS))

  ! Initialise the first section - it is the unnamed one.
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
  type(t_parlist), intent(inout) :: rparlist
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
  type(t_parlist), intent(inout) :: rparlist

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
  type(t_parlist), intent(in) :: rparlist

  ! The section name to look for.
  character(LEN=*), intent(in) :: sname

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
  call sys_toupper (sectionname)

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
  type(t_parlist), intent(inout) :: rparlist

!</inputoutput>

!<input>

  ! The section name to add - without brackets in front and at the end!
  character(LEN=*), intent(in) :: sname

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
  call sys_toupper (sectionname)

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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

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
  call sys_toupper (paramname)

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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

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
  call sys_toupper (paramname)

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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

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
                                           sdefault, isubstring, bdequote)
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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  character(LEN=*), intent(in), optional :: sdefault

  ! OPTIONAL: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  integer, intent(in), optional :: isubstring

  ! OPTIONAL: De-quote the string.
  ! This provides a save way of removing quotation marks around a string in case
  ! the parameter contains exactly one string.
  ! =false: Return the string as it is (standard)
  ! =true: Re-read the string and remove any leading and trailing quotation marks
  !   (if there are any).
  logical, intent(in), optional :: bdequote
!</input>

!<output>

  ! The value of the parameter
  character(LEN=*), intent(out) :: svalue

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
  call sys_toupper (paramname)

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
      call sys_chararraytostring(rsection%p_Rvalues(i)%p_sentry,svalue)
    else
      call sys_chararraytostring(rsection%p_Rvalues(i)%p_SentryList(:,isub),svalue)
    end if
  end if

  if (present(bdequote)) then
    if (bdequote) then
      call sys_dequote(svalue)
    end if
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine parlst_getvalue_string_fetch (rsection, iparameter, svalue,&
                                           bexists, isubstring, bdequote)
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
  type(t_parlstSection), intent(in) :: rsection

  ! The number of the parameter.
  integer, intent(in) :: iparameter

  ! OPTIONAL: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  integer, intent(in), optional :: isubstring

  ! OPTIONAL: De-quote the string.
  ! This provides a save way of removing quotation marks around a string in case
  ! the parameter contains exactly one string.
  ! =false: Return the string as it is (standard)
  ! =true: Re-read the string and remove any leading and trailing quotation marks
  !   (if there are any).
  logical, intent(in), optional :: bdequote
!</input>

!<output>

  ! The value of the parameter
  character(LEN=*), intent(out) :: svalue

  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(out), optional :: bexists

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
    call sys_charArrayToString(rsection%p_Rvalues(iparameter)%p_sentry,svalue)
  else
    call sys_charArrayToString(rsection%p_Rvalues(iparameter)%p_SentryList(:,isub),svalue)
  end if

  if (present(bexists)) bexists = .true.

  if (present(bdequote)) then
    if (bdequote) then
      call sys_dequote(svalue)
    end if
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine parlst_getvalue_string_direct (rparlist, ssectionName, &
                                            sparameter, svalue, sdefault,&
                                            isubstring,bdequote)
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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  character(LEN=*), intent(in), optional :: sdefault

  ! OPTIONAL: The number of the substring to be returned.
  ! =0: returns the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns substring isubstring.
  integer, intent(in), optional :: isubstring

  ! OPTIONAL: De-quote the string.
  ! This provides a save way of removing quotation marks around a string in case
  ! the parameter contains exactly one string.
  ! =false: Return the string as it is (standard)
  ! =true: Re-read the string and remove any leading and trailing quotation marks
  !   (if there are any).
  logical, intent(in), optional :: bdequote

!</input>

!<output>

  ! The value of the parameter
  character(LEN=*), intent(out) :: svalue

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
                                     isubstring,bdequote)

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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  real(SP), intent(in), optional :: fdefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  real(SP), intent(out) :: fvalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  type(t_parlstSection), intent(in) :: rsection

  ! The number of the parameter.
  integer, intent(in) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  real(SP), intent(out) :: fvalue

  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: svalue

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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  real(SP), intent(in), optional :: fdefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  real(SP), intent(out) :: fvalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  real(DP), intent(in), optional :: ddefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  real(DP), intent(out) :: dvalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  type(t_parlstSection), intent(in) :: rsection

  ! The number of the parameter.
  integer, intent(in) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  real(DP), intent(out) :: dvalue

  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: svalue

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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  real(DP), intent(in), optional :: ddefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  real(DP), intent(out) :: dvalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int8_indir (rsection, sparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I8), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I8), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int8_fetch (rsection, iparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The number of the parameter.
  integer, intent(in) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I8), intent(out) :: ivalue

  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: svalue

  svalue = '0'
  call parlst_getvalue_string_fetch (rsection, iparameter, svalue, &
                                     bexists, iarrayindex)
  read(svalue,*) ivalue

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine parlst_getvalue_int8_direct (rparlist, ssectionName, &
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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I8), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I8), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int16_indir (rsection, sparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I16), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I16), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int16_fetch (rsection, iparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The number of the parameter.
  integer, intent(in) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I16), intent(out) :: ivalue

  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: svalue

  svalue = '0'
  call parlst_getvalue_string_fetch (rsection, iparameter, svalue, &
                                     bexists, iarrayindex)
  read(svalue,*) ivalue

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine parlst_getvalue_int16_direct (rparlist, ssectionName, &
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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I16), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I16), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int32_indir (rsection, sparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I32), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I32), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int32_fetch (rsection, iparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The number of the parameter.
  integer, intent(in) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I32), intent(out) :: ivalue

  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: svalue

  svalue = '0'
  call parlst_getvalue_string_fetch (rsection, iparameter, svalue, &
                                     bexists, iarrayindex)
  read(svalue,*) ivalue

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine parlst_getvalue_int32_direct (rparlist, ssectionName, &
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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I32), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I32), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int64_indir (rsection, sparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I64), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I64), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  subroutine parlst_getvalue_int64_fetch (rsection, iparameter, ivalue, &
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
  type(t_parlstSection), intent(in) :: rsection

  ! The number of the parameter.
  integer, intent(in) :: iparameter

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I64), intent(out) :: ivalue

  ! OPTIONAL: Parameter existance check
  ! Is set to TRUE/FALSE, depending on whether the parameter exists.
  logical, intent(out), optional :: bexists

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: svalue

  svalue = '0'
  call parlst_getvalue_string_fetch (rsection, iparameter, svalue, &
                                     bexists, iarrayindex)
  read(svalue,*) ivalue

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine parlst_getvalue_int64_direct (rparlist, ssectionName, &
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
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! OPTIONAL: A default value
  integer(I64), intent(in), optional :: idefault

  ! OPTIONAL: The number of the arrayindex to be returned.
  ! =0: returns the integer directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: returns array index  iarrayindex.
  integer, intent(in), optional :: iarrayindex

!</input>

!<output>

  ! The value of the parameter
  integer(I64), intent(out) :: ivalue

!</output>

!</subroutine>

  ! local variables
  character (LEN=PARLST_LENLINEBUF) :: sdefault,svalue

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
  ! If the parameter exists, it is overwritten.
!</description>

!<inputoutput>

  ! The section where to arr the parameter
  type(t_parlstSection), intent(inout) :: rsection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter
  character(LEN=*), intent(in) :: svalue

  ! OPTIONAL: Number of substrings. This allows a parameter to have
  ! multiple substrings, which can be accessed via the 'isubstring'
  ! parameter in the GET-routines.
  integer, intent(in), optional :: nsubstrings

!</input>

!</subroutine>

    ! local variables
    character(LEN=PARLST_MLNAME) :: paramname
    integer :: i,j

    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call sys_toupper (paramname)

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
      if (associated(rsection%p_Rvalues(i)%p_SentryList)) then
        deallocate(rsection%p_Rvalues(i)%p_SentryList)
        rsection%p_Rvalues(i)%nsize = 0
      end if

      ! Deallocate memory for the value, will be reallocated later.
      if (associated(rsection%p_Rvalues(i)%p_sentry)) then
        deallocate(rsection%p_Rvalues(i)%p_sentry)
      end if

    end if

    rsection%p_Sparameters(i) = paramname
    j = len_trim(svalue)
    allocate(rsection%p_Rvalues(i)%p_sentry(max(1,j)))
    rsection%p_Rvalues(i)%p_sentry(:) = ' '
    call sys_stringtochararray(svalue,rsection%p_Rvalues(i)%p_sentry,j)

    ! Add a list for the substrings if the parameter should have substrings.
    if (present(nsubstrings)) then
      if (nsubstrings .gt. 0) then
        allocate(rsection%p_Rvalues(i)%p_SentryList(max(1,j),nsubstrings))
        rsection%p_Rvalues(i)%p_SentryList(:,:) = ' '
        rsection%p_Rvalues(i)%nsize = nsubstrings
      else
        nullify(rsection%p_Rvalues(i)%p_SentryList)
        rsection%p_Rvalues(i)%nsize = 0
      end if
    else
      nullify(rsection%p_Rvalues(i)%p_SentryList)
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
  type(t_parlist), intent(inout) :: rparlist

!</inputoutput>

!<input>

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The value of the parameter
  character(LEN=*), intent(in) :: svalue

  ! OPTIONAL: Number of substrings. This allows a parameter to have
  ! multiple substrings, which can be accessed via the 'isubstring'
  ! parameter in the GET-routines.
  integer, intent(in), optional :: nsubstrings

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
  type(t_parlstSection), intent(inout) :: rsection

!</inputoutput>

!<input>

  ! The parameter name.
  integer, intent(in) :: iparameter

  ! The new value of the parameter
  character(LEN=*), intent(in) :: svalue

  ! OPTIONAL: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  integer, intent(in), optional :: isubstring

!</input>

!<output>

  ! Optional parameter. Is set to YES/NO, depending on whether
  ! the parameter exists.
  integer, intent(out), optional :: iexists

!</output>

!</subroutine>

  integer :: isub,j

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

  j = len_trim(svalue)
  if ((isub .le. 0) .or. &
      (isub .gt. rsection%p_Rvalues(iparameter)%nsize)) then
    ! Reallocate memory
    deallocate(rsection%p_Rvalues(iparameter)%p_sentry)
    allocate(rsection%p_Rvalues(iparameter)%p_sentry(max(1,j)))
    rsection%p_Rvalues(iparameter)%p_sentry(:) = ' '
    call sys_stringToCharArray(svalue,rsection%p_Rvalues(iparameter)%p_sentry,j)
  else
    ! Check that there is enough memory to save the string.
    if (ubound(rsection%p_Rvalues(iparameter)%p_SentryList,1) .le. j) then
      call parlst_reallocSubVariables(rsection%p_Rvalues(iparameter),j)
    end if
    call sys_stringToCharArray(svalue,rsection%p_Rvalues(iparameter)%p_SentryList(:,isub),j)
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
  type(t_parlstSection), intent(inout) :: rsection

!</inputoutput>

!<input>

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The new value of the parameter
  character(LEN=*), intent(in) :: svalue

  ! OPTIONAL: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  integer, intent(in), optional :: isubstring

!</input>

!</subroutine>

  ! local variables
  integer :: i,isub,j
  character(LEN=PARLST_MLNAME) :: paramname

  ! Create the upper-case parameter name
  paramname = adjustl(sparameter)
  call sys_toupper (paramname)

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

    j = len_trim(svalue)
    if ((isub .le. 0) .or. (isub .gt. rsection%p_Rvalues(i)%nsize)) then
      ! Reallocate memory
      deallocate(rsection%p_Rvalues(i)%p_sentry)
      allocate(rsection%p_Rvalues(i)%p_sentry(max(1,j)))
      rsection%p_Rvalues(i)%p_sentry(:) = ' '
      call sys_stringToCharArray(svalue,rsection%p_Rvalues(i)%p_sentry,j)
    else
      ! Check that there is enough memory to save the string.
      if (ubound(rsection%p_Rvalues(i)%p_SentryList,1) .le. j) then
        call parlst_reallocSubVariables(rsection%p_Rvalues(i),j)
      end if
      call sys_stringToCharArray(svalue,rsection%p_Rvalues(i)%p_SentryList(:,isub),j)
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
  type(t_parlist), intent(inout) :: rparlist

!</inputoutput>

!<input>

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name.
  character(LEN=*), intent(in) :: sparameter

  ! The new value of the parameter
  character(LEN=*), intent(in) :: svalue

  ! OPTIONAL: The number of the substring to be changed.
  ! =0: changes the string directly behind the '=' sign in the line
  !     'name=value'.
  ! >0: changes substring isubstring.
  integer, intent(in), optional :: isubstring

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
  integer, intent(in) :: iunit

  ! The string where to write data to
  character(LEN=*), intent(out) :: sdata

  ! Length of the output
  integer, intent(out) :: ilinelen

  ! Status of the reading process. Set to a value <> 0 if the end
  ! of the file is reached.
  integer, intent(out) :: ios

  ! local variables
  character :: c

  sdata = ''
  ilinelen = 0

  ! Read the data - as long as the line/file does not end.
  do

    ! Read a character.
    ! Unfortunately, Fortran forces me to use this dirty GOTO
    ! to decide processor-independently whether the line or
    ! the record ends.
    read (unit=iunit, fmt='(A1)', iostat=ios, advance='NO',&
          end=10, eor=20) c

    ! Do not do anything in case of an error
    if (ios .eq. 0) then

      ilinelen = ilinelen + 1
      sdata (ilinelen:ilinelen) = c

    end if

    ! Proceed to next character
    cycle

    ! End of file.
10  ios = -1
    exit

    ! End of record = END OF LINE.
20  ios = 0
    exit

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
  character(LEN=*), intent(in) :: sdata

  ! The typ of the line
  integer, intent(out) :: ityp

  ! input: =0: parse line as parameter. isubstring is changed to a value > 0
  !            is the parameter has multiple values attached.
  !        >0: parse line as substring of a multi-valued parameter, not
  !            containing a leading 'name='.
  ! output: If the 'headline' of a multi-valued parameter is read, isubstring is
  !         changed to the number of substrings (the k in 'name(k)=...').
  !         Otherwise unchanged.
  integer, intent(inout) :: isubstring

  ! Line number
  integer, intent(in) :: ilinenum

  ! Section name, if it is a section
  character(LEN=*), intent(inout) :: ssecname

  ! Parameter name, if it is a parameter
  character(LEN=*), intent(inout) :: sparamname

  ! Parameter value, if it is a parameter
  character(LEN=*), intent(inout) :: svalue

  ! OPTIONAL: Filename of the file to be parsed.
  ! Will be printed out in error messages if present.
  character(LEN=*), intent(in), optional :: sfilename

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

        if ((j1 .eq. 0) .or. (j2 .le. j1) .or. (i .le. j1)) then

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
            ! Oh, only one line. User wants to cheat :-)
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

        ! Get the parameter value. Do not get a parameter name; there is none.
        svalue = adjustl(sbuf(1:ltr))

      end if

    end if

  contains

    ! Sub-subroutine: find the length of the line, removing comments
    ! at the end.

    subroutine linelength (sdata, l)

    ! The string to parse. Must not be ''!
    character(LEN=*), intent(in) :: sdata

    ! The index of the last character without any comment at the end.
    integer, intent(out) :: l

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

  subroutine parlst_readfromfile (rparlist, sfile, sdirectory, bexpandVars)

!<description>

  ! This routine parses a text file for data of the INI-file form.
  ! sfile must be the name of a file on the hard disc.
  ! The file may have references to subfiles in the unnamed section:
  ! If there is a parameter "simportdatafiles(n)" present in the
  ! unnamed section, the routine expects a list of filenames in this
  ! parameter. All the files listed there are read in as well.
  ! Parameters at the end of the master file will overwrite
  ! the parameters from the files in simportdatafiles.
  !
  ! Sub-files specified in the main file are searched in the following
  ! directories:
  ! 1.) sdirectory (if specified)
  ! 2.) the directory that contains sfile (if sfile specifies
  !     a directory)
  ! 3.) current directory
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
  type(t_parlist), intent(inout) :: rparlist

!</inputoutput>

!<input>

  ! The filename of the file to read.
  character(LEN=*), intent(in) :: sfile

  ! OPTIONAL: Directory containing other data files in case
  ! the main file sfile contains references to subfiles.
  character(LEN=*), intent(in), optional :: sdirectory

  ! OPTIONAL: Expand references to subvariables.
  ! TRUE: Subvariables like "%{section.varname} are expanded to actual values.
  !       This is the standard setting.
  ! FALSE: Subvariables are left as they are.
  logical, intent(in), optional :: bexpandVars
!</input>

    ! local variables
    integer :: iidxsubfiles,icurrentsubfile
    character(LEN=PARLST_LENLINEBUF), dimension(:), pointer :: p_Ssubfiles,p_SsubfilesTemp
    character(LEN=PARLST_LENLINEBUF) :: sstring,smainfile
    integer :: nsubfiles,nnewsubfiles,j,ilensubf
    logical :: bexists,bmainpath,babsolute
    character(LEN=PARLST_LENLINEBUF) :: smainpath,sfilepath,sfilename

    ! Filename/path of the master dat file.
    ! Search at first in the specified path.
    bexists = .false.
    if (present(sdirectory)) then
      smainfile = trim(sdirectory)//"/"//sfile
      inquire(file=smainfile, exist=bexists)
    end if

    if (.not. bexists) then
      smainfile = sfile
      inquire(file=smainfile, exist=bexists)
    end if

    if (.not. bexists) then
      ! Cancel if the file does not exist.
      return
    end if

    ! Get the main path/name of the file.
    call io_pathExtract (smainfile, smainpath)
    bmainpath = smainpath .ne. ""

    ! Create a list of files to be read.
    ! They contain filenames including the directory.
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

      call io_pathExtract (trim(p_Ssubfiles(icurrentsubfile)), sfilepath, sfilename, babsolute)

      if (babsolute) then
        ! Ok, we have an absolute path given. Test it.
        sstring = trim(p_Ssubfiles(icurrentsubfile))
        inquire(file=sstring, exist=bexists)
      else
        ! Path is relative -- a little bit more complicated.
        ! Is there a directory given?
        if (sfilepath .ne. "") then
          ! Directory specified. We add "sdirectory" if we have it.
          if (present(sdirectory)) then
            sstring = trim(sdirectory)//"/"//trim(sfilepath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if

          if (.not. bexists) then
            ! No, not there. Then take the path directly.
            sstring = trim(sfilepath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if

          if (bmainpath .and. (.not. bexists)) then
            ! No, not there. Add the master directory and test there.
            sstring = trim(smainpath)//"/"//trim(sfilepath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if
        else
          ! No directory given. Then we search in the directory
          ! of the master file...
          if (bmainpath .and. (.not. bexists)) then
            sstring = trim(smainpath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if

          ! And in the current directory.
          if (.not. bexists) then
            sstring = trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if
        end if
      end if

      if (bexists) then
        call parlst_readfromsinglefile (rparlist, sstring, .false., .false.)

        ! Replace the filename with the string including the path.
        ! Then the actual read process at the end of the routine can be
        ! handled easier.
        p_Ssubfiles(icurrentsubfile) = trim(sstring)

      else
        call output_line ('Specified data-subfile does not exist: '//sstring, &
            OU_CLASS_WARNING,OU_MODE_STD,'parlst_readfromfile')
      end if

      ! Check if there is a parameter "simportdatafiles" available in the
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
          call parlst_getvalue_string(rparlist%p_Rsections(1), iidxsubfiles, sstring,bdequote=.true.)
          if (trim(sstring) .ne. "") then
            ! Append the data.
            allocate(p_SsubfilesTemp(nsubfiles+1))
            p_SsubfilesTemp(1:nsubfiles) = p_Ssubfiles(:)
            deallocate(p_Ssubfiles)
            p_Ssubfiles => p_SsubfilesTemp
            nullify(p_SsubfilesTemp)

            ! Expand subvariables and environment variables here.
            ! This point is independent of a parameter bexpandVars
            ! as subfiles may otherwise not be found!
            call parlst_expandEnvVariable(sstring)
            call parlst_expandSubvariable(rparlist,sstring)

            p_Ssubfiles(nsubfiles+1) = sstring
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
                p_Ssubfiles(nsubfiles+j),isubstring=j,bdequote=.true.)

            ! Expand subvariables and environment variables here.
            ! This point is independent of a parameter bexpandVars
            ! as subfiles may otherwise not be found!
            call parlst_expandEnvVariable(p_Ssubfiles(nsubfiles+j))
            call parlst_expandSubvariable(rparlist,p_Ssubfiles(nsubfiles+j))
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
    ilensubf = 1
    do icurrentsubfile = 2,nsubfiles
      sstring = trim(p_Ssubfiles(icurrentsubfile))
      ilensubf = max(ilensubf,len_trim(sstring))
      inquire(file=sstring, exist=bexists)

      if (bexists) then
        ! Read the sub-files...
        ! Do not yet expand variables.
        call parlst_readfromsinglefile (rparlist, sstring, .true., .false.)
      end if
    end do

    ! ... and the master
    icurrentsubfile = 1
    sstring = trim(p_Ssubfiles(icurrentsubfile))
    inquire(file=sstring, exist=bexists)

    if (bexists) then
      ! Do not yet expand variables.
      call parlst_readfromsinglefile (rparlist, sstring, .true., .false.)
    end if

    if (nsubfiles .gt. 1) then
      ! There have a couple of subfiles been read from disc.
      ! All subfiles can be found in p_Ssubfiles.

      ! Incorporate the complete list to the parameter "SIMPORTDATAFILES".
      if (associated(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList)) then
        deallocate(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList)
      end if

      allocate(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList(ilensubf+2,nsubfiles-1))
      rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList(:,:) = ' '
      do icurrentsubfile = 1,nsubfiles-1
        call sys_stringToCharArray('"'//trim(p_Ssubfiles(1+icurrentsubfile))//'"',&
          rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList(:,icurrentsubfile));
      end do
      rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%nsize = nsubfiles-1
    end if

    ! Release memory, finish
    deallocate(p_Ssubfiles)

    ! Now expand all subvariables and environment variables to the actual values.
    if (.not. present(bexpandVars)) then
      call parlst_expandEnvVariables(rparlist)
      call parlst_expandSubvars(rparlist)
    else if (bexpandVars) then
      call parlst_expandEnvVariables(rparlist)
      call parlst_expandSubvars(rparlist)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine parlst_readfromsinglefile (rparlist, sfilename, bimportSections, bexpandVars)

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
  type(t_parlist), intent(inout) :: rparlist

!</inputoutput>

!<input>

  ! The filename of the file to read.
  character(LEN=*), intent(in) :: sfilename

  ! TRUE: Import all sections in the DAT file.
  ! FALSE: Import only the main (unnamed) section and ignore all other
  ! sections.#
  logical, intent(in) :: bimportSections

  ! OPTIONAL: Expand references to subvariables.
  ! TRUE: Subvariables like "%{section.varname} are expanded to actual values.
  !       This is the standard setting.
  ! FALSE: Subvariables are left as they are.
  logical, intent(in), optional :: bexpandVars

!</input>

    ! local variables
    integer :: iunit,ios,isbuflen,ityp,ilinenum,isubstring,nsubstrings,iparpos
    type(t_parlstSection), pointer :: p_currentsection
    character(LEN=PARLST_LENLINEBUF) :: sdata
    character(LEN=PARLST_MLSECTION) :: ssectionname
    character(LEN=PARLST_MLNAME) :: sparname
    character(LEN=PARLST_LENLINEBUF) :: svalue

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

    if (.not. present(bexpandVars)) then
      ! Expand all subvariables and environment variables to the actual values.
      call parlst_expandEnvVariables(rparlist)
      call parlst_expandSubvars(rparlist)
    else if (bexpandVars) then
      ! Expand all subvariables and environment variables to the actual values.
      call parlst_expandEnvVariables(rparlist)
      call parlst_expandSubvars(rparlist)
    end if

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
  type(t_parlist), intent(in) :: rparlist
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
  character(len=PARLST_LENLINEBUF) :: sbuf

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
          call sys_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call appendString(p_sbuf,ilength,&
              trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
              //"="//trim(sbuf))
        else
          ! Write "name(icount)="
          call appendString(p_sbuf,ilength,&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//trim(sys_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          do ientry = 1,icount
            call sys_charArrayToString(&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),sbuf)
            call appendString(p_sbuf,ilength,trim(sbuf))
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
    integer, intent(in) :: nsize

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
    integer, intent(inout) :: iconfigLength

    ! The string to be added.
    character(LEN=*), intent(in) :: sstring

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
  type(t_parlist), intent(in) :: rparlist
!</input>

!</subroutine>


  integer :: isection,ivalue,ientry,icount
  character(len=PARLST_LENLINEBUF) :: sbuf

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
          call sys_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call output_line(&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"="//trim(sbuf))
        else
          ! Write "name(icount)="
          call output_line(&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//trim(sys_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          do ientry = 1,icount
            call sys_charArrayToString(&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
                   p_SentryList(:,ientry),sbuf)
            call output_line(trim(sbuf))
          end do
        end if

      end do ! ivalue

    end do ! isection

  end subroutine

  ! ***************************************************************************

!<function>

  integer function parlst_findvalue_indir (rsection, sparameter, svalue) &
               result (isubstring)

!<description>
  ! Checks whether the parameter sparameter in the section rsection
  ! has the given value svalue.
!</description>

!<result>
  ! The index of the substring in the parameter sparameter which has value
  ! svalue or =-1, if the parameter does not exist within the section.
!</result>

!<input>

  ! The section where to search for the parameter
  type(t_parlstSection), intent(in) :: rsection

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! The value to search for
  character(LEN=*), intent(in) :: svalue

!</input>

!</function>

    ! local variables
    integer :: idx
    character(LEN=PARLST_MLNAME) :: paramname
    character(len(svalue)) :: sbuf
    
    if (sparameter .eq. '') then
      call output_line ('Empty parameter name!', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_lookfor_indir')
      call sys_halt()
    end if
    
    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call sys_toupper (paramname)
    
    ! Get the parameter index into 'idx', finish.
    call parlst_fetchparameter(rsection, paramname, idx)
    
    ! Check if value svalue exists in some substring and return its
    ! index; of the value does not exists return -1
    if (idx .eq. 0) then
      call sys_charArrayToString(&
          rsection%p_Rvalues(idx)%p_sentry, sbuf)
      if (trim(sbuf) .eq. trim(svalue)) then
        isubstring = 0
      else
        isubstring = -1
      end if
    else
      do isubstring = 0, rsection%p_Rvalues(idx)%nsize
        call sys_charArrayToString(&
            rsection%p_Rvalues(idx)%p_SentryList(:,isubstring), sbuf)
        if (trim(sbuf) .eq. trim(svalue)) then
          return
        end if
      end do
      
      ! We did not find the desired value
      isubstring = -1
    end if
    
  end function

  ! ***************************************************************************

!<function>

  integer function parlst_findvalue_direct (rparlist, ssectionName, sparameter, svalue) &
               result (isubstring)

!<description>
  ! Checks whether the parameter sparameter in the section ssectionname
  ! in the parameter list rparlist has the given value.
!</description>

!<result>
  ! The index of the substring in the parameter sparameter which has value
  ! svalue or =-1, if the parameter does not exist within the section.
!</result>

!<input>

  ! The parameter list.
  type(t_parlist), intent(in) :: rparlist

  ! The section name - '' identifies the unnamed section.
  character(LEN=*), intent(in) :: ssectionName

  ! The parameter name to search for.
  character(LEN=*), intent(in) :: sparameter

  ! The value to search for
  character(LEN=*), intent(in) :: svalue

!</input>

!</function>

    ! local variables
    integer :: idx
    type(t_parlstSection), pointer :: p_rsection
    character(len(svalue)) :: sbuf

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call output_line ('Parameter list not initialised!', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_findvalue_direct')
      call sys_halt()
    end if
    
    ! Get the section
    call parlst_querysection(rparlist, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) then
      call output_line ('Section not found: '//trim(ssectionName), &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_findvalue_direct')
      return
    end if

    ! Get the parameter index
    idx = parlst_queryvalue_indir (p_rsection, sparameter)
    
    ! Check if value svalue exists in some substring and return its
    ! index; of the value does not exists return -1
    if (idx .eq. 0) then
      call sys_charArrayToString(&
          p_rsection%p_Rvalues(idx)%p_sentry, sbuf)
      if (trim(sbuf) .eq. trim(svalue)) then
        isubstring = 0
      else
        isubstring = -1
      end if
    else
      do isubstring = 0, p_rsection%p_Rvalues(idx)%nsize
        call sys_charArrayToString(&
            p_rsection%p_Rvalues(idx)%p_SentryList(:,isubstring), sbuf)
        if (trim(sbuf) .eq. trim(svalue)) then
          return
        end if
      end do
      
      ! We did not find the desired value
      isubstring = -1
    end if
    
  end function

  ! ***************************************************************************

!<subroutine>

  subroutine parlst_expandEnvVariables(rparlist)

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
  type(t_parlist), intent(inout) :: rparlist
!</inputoutput>

!</subroutine>

    integer :: isection,ivalue,ientry,icount,j
    character(len=PARLST_LENLINEBUF) :: sbuf

    if (rparlist%isectionCount .eq. 0) then
      call output_line ('Parameter list not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_expandEnvVariables')
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
          call sys_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call parlst_expandEnvVariable(sbuf)

          j = len_trim(sbuf)
          deallocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry)
          allocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(j))
          rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(:) = ' '
          call sys_stringToCharArray(sbuf,&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,j)
        else
          ! Loop through the subvalues.
          do ientry = 1,icount
            ! Expand the value if is refers to subvalues.
            call sys_charArrayToString(&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),sbuf)
            call parlst_expandEnvVariable(sbuf)

            ! Reallocate before writing back if necessary
            j = len_trim(sbuf)
            if (j .gt. ubound(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList,1)) then
              call parlst_reallocSubVariables(rparlist%p_Rsections(isection)%p_Rvalues(ivalue),j)
            end if
            call sys_stringToCharArray(sbuf,&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),j)
          end do
        end if

      end do ! ivalue

    end do ! isection

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine parlst_expandEnvVariable(sbuffer)

  !<description>
    ! This subroutine recursively expands all environment variables in the given
    ! string sbuffer.
  !</description>

  !<input>
    ! string; all environment variables in here are replaced
    character(len=*), intent(inout) :: sbuffer
  !</input>

  !<errors>
    ! ERR_CONV_ENV_VAR_NOT_FOUND (critical)
  !</errors>
!</subroutine>

    ! flag
    logical :: bfoundInEnv

    ! start and end position of variable
    integer(I32) :: istartPos, istopPosRelative

    ! variable to expand environment variable on-the-fly to if found
    character(len=SYS_STRLEN) :: sauxEnv

    ! Buffer for the result
    character(len=PARLST_LENLINEBUF) :: sresult

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

      bfoundInEnv = .false.
      ! Retrieve value of environment variable
      ! (Do not forget to cut the dollar sign.)
      bfoundInEnv = &
          sys_getenv_string(trim(&
              sresult(istartPos + 1 : istartPos + istopPosRelative - 1)), sauxEnv)
      if (bfoundInEnv) then
        ! Replace environment variable by its content
        sresult = sresult(1:istartPos-1) // &
                  trim(sauxEnv) // &
                  trim(sresult(istartPos + istopPosRelative:))
      else
        call output_line ('Environment variable <'//&
            trim(sresult(istartPos + 1 : istartPos + istopPosRelative - 1))//&
            '> not found!',&
            OU_CLASS_ERROR,OU_MODE_STD,'parlst_expandEnvVariables')
        call sys_halt()
      endif

      ! check for next $ character
      istartPos = index(sresult, "$")
    enddo

    ! Replace by the result
    sbuffer = sresult

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine parlst_expandSubvars(rparlist)

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
  type(t_parlist), intent(inout) :: rparlist
!</inputoutput>

!</subroutine>

    integer :: isection,ivalue,ientry,icount,j
    character(len=PARLST_LENLINEBUF) :: sbuf

    if (rparlist%isectionCount .eq. 0) then
      call output_line ('Parameter list not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_expandSubvars')
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
          call sys_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call parlst_expandSubvariable(rparlist,sbuf)

          j = len_trim(sbuf)
          deallocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry)
          allocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(j))
          rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(:) = ' '
          call sys_stringToCharArray(sbuf,&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,j)
        else
          ! Loop through the subvalues.
          do ientry = 1,icount
            ! Expand the value if is refers to subvalues.
            call sys_charArrayToString(&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),sbuf)
            call parlst_expandSubvariable(rparlist,sbuf)

            ! Reallocate before writing back if necessary
            j = len_trim(sbuf)
            if (j .gt. ubound(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList,1)) then
              call parlst_reallocSubVariables(rparlist%p_Rsections(isection)%p_Rvalues(ivalue),j)
            end if
            call sys_stringToCharArray(sbuf,&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),j)
          end do

        end if

      end do ! ivalue

    end do ! isection

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine parlst_findSubvariable(sstring,istart,iend,ssection,sname,ivalue)

!<description>
  ! Searchs in the string sstring for the first occurance of a subvariable.
  ! A subvariable has the form "%{NAME}" or "{SECTION.NAME}"
  ! or "%{NAME:INDEX}" or "%{SECTION.NAME:INDEX}"
  ! depending on whether it is part of the main section or not.
!</description>

!<input>
  ! The string where a subvariable is searched.
  character(len=*), intent(in) :: sstring
!</input>

!<output>
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
!</output>

!</subroutine>

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
    istrlen = len(sstring)
    do i=1,istrlen
      if (sstring(i:i) .eq. "%") then
        ! Did we already found the "%"?
        if (bstartfound) then
          ! That is our escape sequence. Do not do anything, just
          ! return to 'normal' mode.
          bstartfound = .false.
        else
          ! That is probably the beginning of a token.
          bstartfound = .true.
        end if

        ! Go on.
        cycle
      end if

      ! The next things only execute if we are close to a token...
      if (bstartfound) then
        if (sstring(i:I) .eq. "{") then
          ! Yes, that is a token.
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

              ! That is it.
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

  ! ***************************************************************************

!<subroutine>

  subroutine parlst_expandSubvariable(rparlist,sstring)

!<description>
  ! Expands the subvariables in sstring to a fully qualified string.
!</description>

!<input>
  ! The parameter list containing the variables that can be used
  ! as subvariables.
  type(t_parlist), intent(in) :: rparlist
!</input>

!<inputoutput>
  ! The string to be expanded. Receives the expanded string upon return.
  character(len=*), intent(inout) :: sstring
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: istart,iend,ivalue
    character(len=PARLST_MLSECTION) :: ssection
    character(len=PARLST_MLNAME) :: sname
    character(len=len(sstring)) :: sbuffer
    character(len=PARLST_LENLINEBUF) :: sdata

    ! Repeat until we found all subvariables
    istart = 1
    iend = 0
    do while (istart .ne. 0)

      ! Find the first subvariable
      call parlst_findSubvariable(sstring,istart,iend,ssection,sname,ivalue)

      if (istart .ne. 0) then
        ! Copy the string to the buffer
        sbuffer = sstring

        ! Now copy back and replace the variable by the stuff from the
        ! parameter list.
        if (ivalue .eq. 0) then
          call parlst_getvalue_string (rparlist, ssection, sname, sdata, bdequote=.true.)
        else
          call parlst_getvalue_string (rparlist, ssection, sname, sdata, &
              isubstring=ivalue, bdequote=.true.)
        end if
        sstring = sbuffer(1:istart-1)//trim(sdata)//sbuffer(iend+1:)

      end if

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine parlst_dumpToFile(rparlist, sfilename, cflag)

!<description>
   ! This subroutine dumps a given parameter list into a text file in
   ! INI-file form. Note that no additional comments are exported but
   ! just the tuples `parameter = value` in the corresponding sections.
!</description>

!<input>

  ! The parameter list.
  type(t_parlist), intent(in) :: rparlist

  ! The name of the output file
  character(*), intent(in) :: sfilename

  ! mode: SYS_APPEND or SYS_REPLACE
  integer, intent(in) :: cflag

!</input>
!</subroutine>

    ! local variables
    character(len=PARLST_LENLINEBUF) :: sbuf
    integer :: iunit,isection,ivalue,ientry,icount

    if (rparlist%isectionCount .eq. 0) then
      call output_line ('Parameter list not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_dumpToFile')
      call sys_halt()
    end if

    ! Open file for output
    call io_openFileForWriting(sfilename, iunit, cflag, bformatted=.true.)
    if (iunit .eq. -1) then
      call output_line ('Unable to open file for output!', &
          OU_CLASS_ERROR,OU_MODE_STD,'parlst_dumpToFile')
      call sys_halt()
    end if

    ! Loop through all sections
    do isection = 1,rparlist%isectionCount

      ! Append the section name. May be empty for the unnamed section,
      ! which is always the first one.
      if (isection .gt. 1) then
        ! Empty line before
        write(iunit,*)
        write(iunit,*) '['//trim(rparlist%p_Rsections(isection)%ssectionName)//']'
      end if

      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount

        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        if (icount .eq. 0) then
          ! Write "name=value"
          call sys_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          write(iunit,*)&
              trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
              //" = "//trim(sbuf)
        else
          ! Write "name(icount)="
          write(iunit,*)&
              trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
              //"("//trim(sys_siL(icount, 10))//") = "
          ! Write all the entries of that value, one each line.
          do ientry = 1,icount
            call sys_charArrayToString(&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
                p_SentryList(:,ientry),sbuf)
            write(iunit,*) trim(sbuf)
          end do
        end if

      end do ! ivalue

    end do ! isection

    ! Close file
    close(iunit)

  end subroutine parlst_dumpToFile

end module
