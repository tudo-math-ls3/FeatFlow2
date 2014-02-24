!##############################################################################
!# ****************************************************************************
!# <name> convergencetable </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a convergence table that can be used to
!# generate pretty convergence tables in Latex
!#
!# The following routines are provided by this module:
!#
!# 1.) ctab_init
!#     -> Initialises an empty convergence table
!#
!# 2.) ctab_addValue
!#     -> Adds a new value to column of the convergence table
!#
!# 3.) ctab_setPrecision
!#     -> Sets the precision for a column for tex output
!#
!# 4.) ctab_setScientific
!#     -> Sets scientific notation for a column for tex output
!#
!# 5.) ctab_setTexCaption
!#     -> Sets the tex caption of a single column for tex output
!#
!# 6.) ctab_setTexTableCaption
!#     -> Sets the tex caption of an entire table for tex output
!#
!# 7.) ctab_setTexTableLabel
!#     -> Sets the tex label of an entire table for tex output
!#
!# 8.) ctab_setTexFormat
!#     -> Sets the tex format of a single column for tex output
!#
!# 9.) ctab_getColumn
!#     -> Get pointer to column
!#
!# 10.) ctab_getCell
!#      -> Get pointer to cell
!#
!# 11.) ctab_swapColumns
!#      -> Swaps two columns
!#
!# 12.) ctab_outputText
!#      -> Writes table to screen in human-readable format
!#
!# 13.) ctab_outputTex
!#      -> Writes table to screen in tex format
!#
!# 14.) ctab_evalConvergenceRate
!#      -> Evalualtes the convergence rates for a data column
!#
!# 15.) ctab_done
!#      -> Releases the convergence table
!# </purpose>
!##############################################################################

module convergencetable

!$use omp_lib
  use fsystem
  use genoutput
  use io
  use storage

  implicit none

  private

  public :: t_convergenceTable
  public :: ctab_init
  public :: ctab_addValue
  public :: ctab_setPrecision
  public :: ctab_setScientific
  public :: ctab_setTexCaption
  public :: ctab_setTexTableCaption
  public :: ctab_setTexTableLabel
  public :: ctab_setTexFormat
  public :: ctab_getColumn
  public :: ctab_getCell
  public :: ctab_swapColumns
  public :: ctab_outputText
  public :: ctab_outputTex
  public :: ctab_evalConvergenceRate
  public :: ctab_done

!<constants>

!<constantblock description="Constants for convergence rate evaluation">

  ! Compute reduction rate
  integer, parameter, public :: CTAB_REDUCTION_RATE      = 0

  ! Compute logarithm of reduction rate
  integer, parameter, public :: CTAB_REDUCTION_RATE_LOG2 = 1

!</constantblock>

!</constants>

!<types>

!<typeblock>
  ! Structure defining a cell of a convergence table
  type t_cell
    
    private

    ! Pointer to previous and next cell
    type(t_cell), pointer :: p_rprevCell => null()
    type(t_cell), pointer :: p_rnextCell => null()
    
    ! Content of cell
    character(LEN=8) :: svalue = ''

    ! Type of content
    integer :: ctype = ST_NOHANDLE

  end type t_cell
!</typeblock>

!<typeblock>
  ! Structure defining a column of a convergence table
  type t_column
    
    private

    ! Column name
    character(LEN=SYS_STRLEN) :: skey = ''

    ! Number of decimal digits
    integer :: idigits = 8

    ! Scientific flag
    logical :: bscientific = .false.

    ! Tex caption string
    character(LEN=SYS_STRLEN) :: stexCaption = ''

    ! Tex format string
    character(LEN=SYS_STRLEN) :: stexFormat = 'l'

    ! Pointer to previous and next column
    type(t_column), pointer :: p_rprevColumn => null()
    type(t_column), pointer :: p_rnextColumn => null()

    ! Pointer to first and last cell in column
    type(t_cell), pointer :: p_rfirstCell => null()
    type(t_cell), pointer :: p_rlastCell => null()

  end type t_column
!</typeblock>

!<typeblock>

  ! Structure defining a convergence table
  type t_convergenceTable

    private

    ! Tex caption string
    character(LEN=SYS_STRLEN) :: stexCaption = ''

    ! Tex label string
    character(LEN=SYS_STRLEN) :: stexLabel = ''
    
    ! Pointer to first and last column in table
    type(t_column), pointer :: p_rfirstColumn => null()
    type(t_column), pointer :: p_rlastColumn => null()

  end type t_convergenceTable
!</typeblock>

!</types>

  interface ctab_addValue
    module procedure ctab_addValue_single
    module procedure ctab_addValue_double
    module procedure ctab_addValue_int8
    module procedure ctab_addValue_int16
    module procedure ctab_addValue_int32
    module procedure ctab_addValue_int64
  end interface ctab_addValue

  interface ctab_getColumn
    module procedure ctab_getColumnByKey
    module procedure ctab_getColumnByNumber
  end interface ctab_getColumn

  interface ctab_swapColumns
    module procedure ctab_swapColumnsByKey
    module procedure ctab_swapColumnsByNumber
    module procedure ctab_swapColumnsByObject
  end interface ctab_swapColumns

  interface ctab_evalConvergenceRate
    module procedure ctab_evalConvRateByKey
    module procedure ctab_evalConvRateByNumber
    module procedure ctab_evalConvRateByObject
  end interface ctab_evalConvergenceRate

contains

  !************************************************************************

!<subroutine>

  subroutine ctab_init(rtable)

!<description>
    ! This subroutine initialises an empty convergence table
!</description>

!<output>
    ! A convergence table
    type(t_convergenceTable), intent(out) :: rtable
!</output>

!</subroutine>

    ! Nullify pointers
    nullify(rtable%p_rfirstColumn)
    nullify(rtable%p_rlastColumn)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_addValue_single(rtable,skey,value,p_result)

!<description>
    ! This subroutine adds a single precision value to the column with
    ! key SKEY. If no column with key SKEY exists then a new column is
    ! created at the end of the table.
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Single precision value
    real(SP), intent(in) :: value
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!<output>
    ! OPTIONAL: Pointer to the new cell
    type(t_cell), pointer, optional :: p_result
!</output>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell

    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)

    ! Create new column if it does not exist
    if (.not.associated(p_rcolumn)) then
      allocate(p_rcolumn)
      p_rcolumn%skey=skey
      call ctab_addColumn(rtable,p_rcolumn)
    end if

    ! Add new cell and transfer value
    allocate(p_rcell)
    p_rcell%svalue = transfer(value,p_rcell%svalue)
    p_rcell%ctype = ST_SINGLE

    call ctab_addCell(p_rcolumn,p_rcell)
    if (present(p_result)) p_result => p_rcell

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_addValue_double(rtable,skey,value,p_result)

!<description>
    ! This subroutine adds a double precision value to the column with
    ! key SKEY. If no column with key SKEY exists then a new column is
    ! created at the end of the table.
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Double precision value
    real(DP), intent(in) :: value
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!<output>
    ! OPTIONAL: Pointer to the new cell
    type(t_cell), pointer, optional :: p_result
!</output>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell

    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)

    ! Create new column if it does not exist
    if (.not.associated(p_rcolumn)) then
      allocate(p_rcolumn)
      p_rcolumn%skey=skey
      call ctab_addColumn(rtable,p_rcolumn)
    end if

    ! Add new cell and transfer value
    allocate(p_rcell)
    p_rcell%svalue = transfer(value,p_rcell%svalue)
    p_rcell%ctype = ST_DOUBLE

    call ctab_addCell(p_rcolumn,p_rcell)
    if (present(p_result)) p_result => p_rcell

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_addValue_int8(rtable,skey,value,p_result)

!<description>
    ! This subroutine adds an integer value to the column with key
    ! SKEY. If no column with key SKEY exists then a new column is
    ! created at the end of the table.
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Integer value
    integer(I8), intent(in) :: value
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!<output>
    ! OPTIONAL: Pointer to the new cell
    type(t_cell), pointer, optional :: p_result
!</output>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell

    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)

    ! Create new column if it does not exist
    if (.not.associated(p_rcolumn)) then
      allocate(p_rcolumn)
      p_rcolumn%skey=skey
      call ctab_addColumn(rtable,p_rcolumn)
    end if

    ! Add new cell and transfer value
    allocate(p_rcell)
    p_rcell%svalue = transfer(value,p_rcell%svalue)
    p_rcell%ctype = ST_INT8

    call ctab_addCell(p_rcolumn,p_rcell)   
    if (present(p_result)) p_result => p_rcell

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_addValue_int16(rtable,skey,value,p_result)

!<description>
    ! This subroutine adds an integer value to the column with key
    ! SKEY. If no column with key SKEY exists then a new column is
    ! created at the end of the table.
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Integer value
    integer(I16), intent(in) :: value
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!<output>
    ! OPTIONAL: Pointer to the new cell
    type(t_cell), pointer, optional :: p_result
!</output>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell

    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)

    ! Create new column if it does not exist
    if (.not.associated(p_rcolumn)) then
      allocate(p_rcolumn)
      p_rcolumn%skey=skey
      call ctab_addColumn(rtable,p_rcolumn)
    end if

    ! Add new cell and transfer value
    allocate(p_rcell)
    p_rcell%svalue = transfer(value,p_rcell%svalue)
    p_rcell%ctype = ST_INT16

    call ctab_addCell(p_rcolumn,p_rcell)
    if (present(p_result)) p_result => p_rcell

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_addValue_int32(rtable,skey,value,p_result)

!<description>
    ! This subroutine adds an integer value to the column with key
    ! SKEY. If no column with key SKEY exists then a new column is
    ! created at the end of the table.
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Integer value
    integer(I32), intent(in) :: value
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!<output>
    ! OPTIONAL: Pointer to the new cell
    type(t_cell), pointer, optional :: p_result
!</output>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell

    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)

    ! Create new column if it does not exist
    if (.not.associated(p_rcolumn)) then
      allocate(p_rcolumn)
      p_rcolumn%skey=skey
      call ctab_addColumn(rtable,p_rcolumn)
    end if

    ! Add new cell and transfer value
    allocate(p_rcell)
    p_rcell%svalue = transfer(value,p_rcell%svalue)
    p_rcell%ctype = ST_INT32

    call ctab_addCell(p_rcolumn,p_rcell)
    if (present(p_result)) p_result => p_rcell

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_addValue_int64(rtable,skey,value,p_result)

!<description>
    ! This subroutine adds an integer value to the column with key
    ! SKEY. If no column with key SKEY exists then a new column is
    ! created at the end of the table.
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Integer value
    integer(I64), intent(in) :: value
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!<output>
    ! OPTIONAL: Pointer to the new cell
    type(t_cell), pointer, optional :: p_result
!</output>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell

    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)

    ! Create new column if it does not exist
    if (.not.associated(p_rcolumn)) then
      allocate(p_rcolumn)
      p_rcolumn%skey=skey
      call ctab_addColumn(rtable,p_rcolumn)
    end if

    ! Add new cell and transfer value
    allocate(p_rcell)
    p_rcell%svalue = transfer(value,p_rcell%svalue)
    p_rcell%ctype = ST_INT64

    call ctab_addCell(p_rcolumn,p_rcell)   
    if (present(p_result)) p_result => p_rcell

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_setPrecision(rtable,skey,idigits)

!<description>
    ! This subroutine sets the precision of the column with key SKEY
    ! to be used for writing e.g. double or single values
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Precision
    integer, intent(in) :: idigits
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>
    
    ! local variables
    type(t_column), pointer :: p_rcolumn
    
    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)
    
    if (associated(p_rcolumn)) then
      p_rcolumn%idigits = idigits
    else
      ! Throw error if column does not exist
      call output_line("Column does not exist!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_setPrecision")
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_setScientific(rtable,skey,bscientific)

!<description>
    ! This subroutine sets the scientific flag of the column with key SKEY
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Scientific flag
    logical, intent(in) :: bscientific
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>
    
    ! local variables
    type(t_column), pointer :: p_rcolumn
    
    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)
    
    if (associated(p_rcolumn)) then
      p_rcolumn%bscientific = bscientific
    else
      ! Throw error if column does not exist
      call output_line("Column does not exist!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_setScientific")
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_setTexCaption(rtable,skey,sstring)

!<description>
    ! This subroutine sets the caption of the column with key SKEY to
    ! be used for tex output
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Tex column caption string
    character(len=*), intent(in) :: sstring
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>
    
    ! local variables
    type(t_column), pointer :: p_rcolumn
    
    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)
    
    if (associated(p_rcolumn)) then
      p_rcolumn%stexCaption = sstring
    else
      ! Throw error if column does not exist
      call output_line("Column does not exist!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_setTexCaption")
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_setTexTableCaption(rtable,sstring)

!<description>
    ! This subroutine sets the caption of the table to be used for tex
    ! output
!</description>

!<input>
    ! Tex table caption string
    character(len=*), intent(in) :: sstring
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    rtable%sTexCaption = sstring

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_setTexTableLabel(rtable,sstring)

!<description>
    ! This subroutine sets the label of the table to be used for tex
    ! output
!</description>

!<input>
    ! Tex table label string
    character(len=*), intent(in) :: sstring
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    rtable%sTexLabel = sstring

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_setTexFormat(rtable,skey,sstring)

!<description>
    ! This subroutine sets the format of the column with key SKEY to
    ! be used for tex output
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Tex column format string
    character(len=*), intent(in) :: sstring
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>
    
    ! local variables
    type(t_column), pointer :: p_rcolumn
    
    ! Get column
    p_rcolumn => ctab_getColumn(rtable,skey)
    
    if (associated(p_rcolumn)) then
      p_rcolumn%stexFormat = sstring
    else
      ! Throw error if column does not exist
      call output_line("Column does not exist!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_setTexFormat")
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<function>

  function ctab_getColumnByKey(rtable,skey) result(p_rcolumn)

!<description>
    ! This functions returns a pointer to the column with key SKEY.
    ! If no column with key SKEY exists then NULL is returned.
!</description>

!<input>
    ! Convergence table
    type(t_convergenceTable), intent(in) :: rtable

    ! Key of the column
    character(len=*), intent(in) :: skey
!</input>

!<result>
    type(t_column), pointer :: p_rcolumn
!</result>

!</function>

    p_rcolumn => rtable%p_rfirstColumn
    do while (associated(p_rcolumn))
      if (trim(p_rcolumn%skey) .eq. skey) return
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do

  end function

  !************************************************************************

!<function>

  function ctab_getColumnByNumber(rtable,inumber) result(p_rcolumn)

!<description>
    ! This functions returns a pointer to the INUMBER-th column.
    ! If no INUMBER-th column exists then NULL is returned.
!</description>

!<input>
    ! Convergence table
    type(t_convergenceTable), intent(in) :: rtable

    ! Number of the column
    integer, intent(in) :: inumber
!</input>

!<result>
    type(t_column), pointer :: p_rcolumn
!</result>

!</function>

    ! local variables
    integer :: i

    i=1
    p_rcolumn => rtable%p_rfirstColumn
    do while (associated(p_rcolumn))
      if (i .eq. inumber) return
      i=i+1
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do

  end function

  !************************************************************************

!<function>

  function ctab_getCell(rcolumn,inumber) result(p_rcell)

!<description>
    ! This functions returns a pointer to the INUMBER-th cell.
    ! If no INUMBER-th exists cell then NULL is returned.
!</description>

!<input>
    ! Column
    type(t_column), intent(in) :: rcolumn

    ! Number of the column
    integer, intent(in) :: inumber
!</input>

!<result>
    type(t_cell), pointer :: p_rcell
!</result>

!</function>

    ! local variables
    integer :: i

    i=1
    p_rcell => rcolumn%p_rfirstCell
    do while (associated(p_rcell))
      if (i .eq. inumber) return
      i=i+1
      p_rcell => p_rcell%p_rnextCell
    end do

  end function

  !************************************************************************

!<subroutine>

  subroutine ctab_swapColumnsByKey(rtable,skey1,skey2)

!<description>
    ! This subroutine swaps the two columns ICOLUMN and JCOLUMN
!</description>

!<input>
    ! Keys of the columns to be swapped
    character(len=*), intent(in) :: skey1,skey2
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn1,p_rcolumn2
    type(t_column), pointer :: p_rprevColumn,p_rnextColumn

    ! Check if both columns are the same
    if (skey2 .eq. skey2) return

    ! Get columns
    p_rcolumn1 => ctab_getColumn(rtable,skey1)
    p_rcolumn2 => ctab_getColumn(rtable,skey2)

    call ctab_swapColumnsByObject(rtable,p_rcolumn1,p_rcolumn2)
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_swapColumnsByNumber(rtable,inumber1,inumber2)

!<description>
    ! This subroutine swaps the two columns INUMBER1 and INUMBER2
!</description>

!<input>
    ! Column numbers to be swapped
    integer, intent(in) :: inumber1,inumber2
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn1,p_rcolumn2
    type(t_column), pointer :: p_rprevColumn1,p_rnextColumn1
    type(t_column), pointer :: p_rprevColumn2,p_rnextColumn2

    ! Check if both columns are the same
    if (inumber1 .eq. inumber2) return

    ! Get columns
    p_rcolumn1 => ctab_getColumn(rtable,inumber1)
    p_rcolumn2 => ctab_getColumn(rtable,inumber2)

    call ctab_swapColumnsByObject(rtable,p_rcolumn1,p_rcolumn2)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_swapColumnsByObject(rtable,p_rcolumn1,p_rcolumn2)

!<description>
    ! This subroutine swaps the two columns p_rcolumns1 and p_rcolumns2
!</description>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable

    ! Columns to be swapped
    type(t_column), pointer :: p_rcolumn1,p_rcolumn2
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rprevColumn1,p_rnextColumn1
    type(t_column), pointer :: p_rprevColumn2,p_rnextColumn2

    ! Check if both columns are the same
    if (associated(p_rcolumn1,p_rcolumn2)) return

    if (associated(p_rcolumn1) .and. associated(p_rcolumn2)) then

      ! Check if either of the columns is first of the table
      if (associated(rtable%p_rfirstColumn,p_rcolumn1)) then
        rtable%p_rfirstColumn => p_rcolumn2
      elseif (associated(rtable%p_rfirstColumn,p_rcolumn2)) then
        rtable%p_rfirstColumn => p_rcolumn1
      end if

      ! Check if either of the columns is first of the table
      if (associated(rtable%p_rlastColumn,p_rcolumn1)) then
        rtable%p_rlastColumn => p_rcolumn2
      elseif (associated(rtable%p_rlastColumn,p_rcolumn2)) then
        rtable%p_rlastColumn => p_rcolumn1
      end if

      ! Save pointers to prev/next columns
      p_rnextColumn1 => p_rcolumn1%p_rnextColumn
      p_rnextColumn2 => p_rcolumn2%p_rnextColumn
      p_rprevColumn1 => p_rcolumn1%p_rprevColumn
      p_rprevColumn2 => p_rcolumn2%p_rprevColumn

      ! Check if both column follow each other
      if (associated(p_rcolumn1,p_rprevColumn2) .and.&
          associated(p_rcolumn2,p_rnextColumn1)) then
        
        ! Adjust prev/next of columns
        if (associated(p_rprevColumn1)) p_rprevColumn1%p_rnextColumn => p_rcolumn2
        p_rcolumn2%p_rprevColumn => p_rprevColumn1
        
        if (associated(p_rnextColumn2)) p_rnextColumn2%p_rprevColumn => p_rcolumn1
        p_rcolumn1%p_rnextColumn => p_rnextColumn2
        
        p_rcolumn1%p_rprevColumn => p_rcolumn2
        p_rcolumn2%p_rnextColumn => p_rcolumn1

      elseif (associated(p_rcolumn2,p_rprevColumn1) .and.&
              associated(p_rcolumn1,p_rnextColumn2)) then

        ! Adjust prev/next of columns
        if (associated(p_rprevColumn2)) p_rprevColumn2%p_rnextColumn => p_rcolumn1
        p_rcolumn1%p_rprevColumn => p_rprevColumn2
        
        if (associated(p_rnextColumn1)) p_rnextColumn1%p_rprevColumn => p_rcolumn2
        p_rcolumn2%p_rnextColumn => p_rnextColumn1
        
        p_rcolumn2%p_rprevColumn => p_rcolumn1
        p_rcolumn1%p_rnextColumn => p_rcolumn2

      else

        ! Adjust prev/next of columns
        if (associated(p_rprevColumn1)) p_rprevColumn1%p_rnextColumn => p_rcolumn2
        if (associated(p_rprevColumn2)) p_rprevColumn2%p_rnextColumn => p_rcolumn1
        p_rcolumn1%p_rprevColumn => p_rprevColumn2
        p_rcolumn2%p_rprevColumn => p_rprevColumn1
        
        if (associated(p_rnextColumn1)) p_rnextColumn1%p_rprevColumn => p_rcolumn2
        if (associated(p_rnextColumn2)) p_rnextColumn2%p_rprevColumn => p_rcolumn1
        p_rcolumn1%p_rnextColumn => p_rnextColumn2
        p_rcolumn2%p_rnextColumn => p_rnextColumn1

      end if

    else 
      ! Throw error if column does not exist
      call output_line("Column does not exist!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_swapColumnsByNumber")
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_outputText(rtable,ioutputUnit)

!<description>
    ! This subroutine writes a convergence table in human-readable
    ! format.  If the optional parameter IUNIT is given, then the
    ! output is redirected to unit IUNIT. Otherwise the standard
    ! output unit is used.
!</description>

!<input>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable

    ! OPTIONAL: output unit
    integer, intent(in), optional :: ioutputUnit
!</input>

!</subroutine>

    ! local variable
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell
    character(len=1024) :: sstring
    integer, dimension(:), allocatable :: IcolumnWidth,IcolumnTab
    integer :: iunit,icolumn,irow,ncolumns,nrows
    real(SP) :: fvalue
    real(DP) :: dvalue
    integer(I8)  :: i8value
    integer(I16) :: i16value
    integer(I32) :: i32value
    integer(I64) :: i64value

    ! Check if explicit output unit is given
    if (present(ioutputUnit)) then
      iunit = ioutputUnit
    else
      iunit = OU_TERMINAL
    end if

    ! Compute the number of columns/rows
    ncolumns = ctab_getNColumns(rtable)
    nrows = ctab_getNRows(rtable)

    ! Compute the column width
    allocate(IcolumnWidth(ncolumns),IcolumnTab(ncolumns+1))
    
    icolumn = 1
    p_rcolumn => rtable%p_rfirstColumn
    do while(associated(p_rcolumn))
      IcolumnWidth(icolumn) = ctab_getColumnWidth(p_rcolumn)
      icolumn = icolumn+1
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do

    ! Compute column tab-stops
    IcolumnTab(1) = 1
    do icolumn=2,ncolumns+1
      IcolumnTab(icolumn) = IcolumnTab(icolumn-1)+IcolumnWidth(icolumn-1)+1
    end do

    ! Write captions
    sstring = ''
    icolumn = 1
    p_rcolumn => rtable%p_rfirstColumn
    do while(associated(p_rcolumn))
      sstring(IcolumnTab(icolumn):) = trim(adjustl(p_rcolumn%skey))
      icolumn = icolumn+1
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do
    write(iunit,fmt='(A)') trim(sstring)

    ! Write content of cells
    do irow=1,nrows
      sstring = ''
      icolumn = 1
      p_rcolumn => rtable%p_rfirstColumn
      do while(associated(p_rcolumn))

        p_rcell => ctab_getCell(p_rcolumn,irow)
        if (associated(p_rcell)) then
          
          ! What type of data are we?
          select case(p_rcell%ctype)
          case (ST_SINGLE)
            fvalue = transfer(p_rcell%svalue,fvalue)
            if (p_rcolumn%bscientific) then
              sstring(IcolumnTab(icolumn):) =&
                  sys_adjustr(sys_sdE(real(fvalue,DP),p_rcolumn%idigits),&
                              IcolumnWidth(icolumn))
            else
              sstring(IcolumnTab(icolumn):) =&
                  sys_adjustr(sys_sd(real(fvalue,DP),p_rcolumn%idigits),&
                              IcolumnWidth(icolumn))
            end if
        
          case (ST_DOUBLE)
            dvalue = transfer(p_rcell%svalue,dvalue)
            if (p_rcolumn%bscientific) then
              sstring(IcolumnTab(icolumn):) =&
                  sys_adjustr(sys_sdE(dvalue,p_rcolumn%idigits),&
                              IcolumnWidth(icolumn))
            else
              sstring(IcolumnTab(icolumn):) =&
                  sys_adjustr(sys_sd(dvalue,p_rcolumn%idigits),&
                              IcolumnWidth(icolumn))

            end if

          case (ST_INT8)
            i8value = transfer(p_rcell%svalue,i8value)
            sstring(IcolumnTab(icolumn):) =&
                sys_adjustr(sys_si(int(i8value),IcolumnWidth(icolumn)),&
                            IcolumnWidth(icolumn))
            
          case (ST_INT16)
            i16value = transfer(p_rcell%svalue,i16value)
            sstring(IcolumnTab(icolumn):) =&
                sys_adjustr(sys_si(int(i16value),IcolumnWidth(icolumn)),&
                            IcolumnWidth(icolumn))
            
          case (ST_INT32)
            i32value = transfer(p_rcell%svalue,i32value)
            sstring(IcolumnTab(icolumn):) =&
                sys_adjustr(sys_si(int(i32value),IcolumnWidth(icolumn)),&
                            IcolumnWidth(icolumn))
                        
          case (ST_INT64)
            i64value = transfer(p_rcell%svalue,i64value)
            sstring(IcolumnTab(icolumn):) =&
                sys_adjustr(sys_sli(i64value,IcolumnWidth(icolumn)),&
                            IcolumnWidth(icolumn))
            
          case default
            ! Throw error if data type is not supported
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"ctab_outputText")
            call sys_halt()
          end select
        end if
        
        icolumn = icolumn+1
        p_rcolumn => p_rcolumn%p_rnextColumn
      end do
      write(iunit,fmt='(A)') trim(sstring)
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_outputTex(rtable,soutputFile)

!<description>
    ! This subroutine writes a convergence table in Tex format.
!</description>

!<input>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable

    ! Name of the Tex file
    character(len=*), intent(in) :: soutputFile
!</input>

!</subroutine>

    ! local variable
    type(t_column), pointer :: p_rcolumn
    type(t_cell), pointer :: p_rcell
    character(len=1024) :: sstring
    integer :: iunit,icolumn,irow,ncolumns,nrows
    real(SP) :: fvalue
    real(DP) :: dvalue
    integer(I8)  :: i8value
    integer(I16) :: i16value
    integer(I32) :: i32value
    integer(I64) :: i64value

    ! Open file for writing
    call io_openFileForWriting(soutputFile,iunit,SYS_REPLACE)

    ! Compute the number of columns/rows
    ncolumns = ctab_getNColumns(rtable)
    nrows = ctab_getNRows(rtable)

    ! Write Tex header
    write(iunit,fmt='(A)') '\begin{table}'
    sstring = '\begin{tabular}[c]{|'
    p_rcolumn => rtable%p_rfirstColumn
    do while(associated(p_rcolumn))
      sstring = trim(sstring)//trim(p_rcolumn%stexFormat)//'|'
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do
    sstring = trim(sstring)//'}'
    write(iunit,fmt='(A)') trim(sstring)
    write(iunit,fmt='(A)') '\hline'
    
    ! Write captions
    sstring = ''
    icolumn = 1
    p_rcolumn => rtable%p_rfirstColumn
    do while(associated(p_rcolumn))
      if (icolumn .lt. ncolumns) then
        sstring = trim(sstring)//trim(adjustl(p_rcolumn%sTexCaption))//'&'
      else
        sstring = trim(sstring)//trim(adjustl(p_rcolumn%sTexCaption))//'\\'
      end if
      icolumn = icolumn+1
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do
    write(iunit,fmt='(A)') trim(sstring)
    write(iunit,fmt='(A)') '\hline'

    ! Write content of cells
    do irow=1,nrows
      sstring = ''
      icolumn = 1
      p_rcolumn => rtable%p_rfirstColumn
      do while(associated(p_rcolumn))

        p_rcell => ctab_getCell(p_rcolumn,irow)
        if (associated(p_rcell)) then
          
          ! What type of data are we?
          select case(p_rcell%ctype)
          case (ST_SINGLE)
            fvalue = transfer(p_rcell%svalue,fvalue)
            if (p_rcolumn%bscientific) then
              sstring = trim(sstring)//&
                  trim(sys_sdE(real(fvalue,DP),p_rcolumn%idigits))
            else
              sstring = trim(sstring)//&
                  trim(sys_sd(real(fvalue,DP),p_rcolumn%idigits))
            end if
        
          case (ST_DOUBLE)
            dvalue = transfer(p_rcell%svalue,dvalue)
            if (p_rcolumn%bscientific) then
              sstring = trim(sstring)//&
                  trim(sys_sdE(dvalue,p_rcolumn%idigits))
            else
              sstring = trim(sstring)//&
                  trim(sys_sd(dvalue,p_rcolumn%idigits))
            end if

          case (ST_INT8)
            i8value = transfer(p_rcell%svalue,i8value)
            sstring = trim(sstring)//trim(sys_si(int(i8value),16))
            
          case (ST_INT16)
            i16value = transfer(p_rcell%svalue,i16value)
            sstring = trim(sstring)//trim(sys_si(int(i16value),16))

          case (ST_INT32)
            i32value = transfer(p_rcell%svalue,i32value)
            sstring = trim(sstring)//trim(sys_si(int(i32value),16))
                       
          case (ST_INT64)
            i64value = transfer(p_rcell%svalue,i64value)
            sstring = trim(sstring)//trim(sys_sli(i64value,32))
            
          case default
            ! Throw error if data type is not supported
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"ctab_outputText")
            call sys_halt()
          end select
        end if
        
        if (icolumn .lt. ncolumns) then
          sstring = trim(sstring)//'&'
        else
          sstring = trim(sstring)//'\\'
        end if
        icolumn = icolumn+1
        p_rcolumn => p_rcolumn%p_rnextColumn
      end do
      write(iunit,fmt='(A)') trim(sstring)
      write(iunit,fmt='(A)') '\hline'
    end do
      
    ! Write Tex footer (caption,label)
    write(iunit,fmt='(A)') '\end{tabular}'
    if (rtable%sTexCaption .ne. '') write(iunit,fmt='(A)')&
        '\caption{'//trim(rtable%sTexCaption)//'}'
    if (rtable%sTexLabel .ne. '') write(iunit,fmt='(A)')&
        '\label{'//trim(rtable%sTexLabel)//'}'
    write(iunit,fmt='(A)') '\end{table}'
    
    ! Close file
    close(iunit)
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_evalConvRateByKey(rtable,skey,cevalType)

!<description>
    ! This subroutine computes the convergence rate for the column
    ! with key SKEY and inserts a new column right after it
!</description>

!<input>
    ! Key of the column
    character(len=*), intent(in) :: skey

    ! Type of convergence rate evaluation. One of the CTAB_xxx constants.
    ! CTAB_REDUCTION_RATE:      Reduction rate
    ! CTAB_REDUCTION_RATE_LOG2: Logarithm of reduction rate
    integer, intent(in) :: cevalType
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>
    
    ! local variables
    type(t_column), pointer :: p_rcolumnData,p_rcolumnConvRate

    ! Get column
    p_rcolumnData => ctab_getColumn(rtable,skey)

    if (associated(p_rcolumnData)) then

      ! Create new column with convergence rates
      allocate(p_rcolumnConvRate)
      p_rcolumnConvRate%skey = trim(skey)//'-convrate'
      call ctab_addColumn(rtable,p_rcolumnConvRate)
      call ctab_evalConvergenceRate(p_rcolumnData,p_rcolumnConvRate,cevalType)
      
      ! Swap columns
      p_rcolumnData => p_rcolumnData%p_rnextColumn
      call ctab_swapColumns(rtable,p_rcolumnData,p_rcolumnConvRate)
    else
      ! Throw error if column does not exist
      call output_line("Column does not exist!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_evalConvRateByKey")
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_evalConvRateByNumber(rtable,inumber,cevalType)

!<description>
    ! This subroutine computes the convergence rate for the INUMBER-th
    ! column and inserts a new column right after it
!</description>

!<input>
    ! Column number
    integer, intent(in) :: inumber

    ! Type of convergence rate evaluation. One of the CTAB_xxx constants.
    ! CTAB_REDUCTION_RATE:      Reduction rate
    ! CTAB_REDUCTION_RATE_LOG2: Logarithm of reduction rate
    integer, intent(in) :: cevalType
!</input>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>
    
    ! local variables
    type(t_column), pointer :: p_rcolumnData,p_rcolumnConvRate

    ! Get column
    p_rcolumnData => ctab_getColumn(rtable,inumber)

    if (associated(p_rcolumnData)) then
      
      ! Create new column with convergence rates
      allocate(p_rcolumnConvRate)
      p_rcolumnConvRate%skey = trim(p_rcolumnData%skey)//'-convrate'
      call ctab_addColumn(rtable,p_rcolumnConvRate)
      call ctab_evalConvergenceRate(p_rcolumnData,p_rcolumnConvRate,cevalType)

      ! Swap columns
      p_rcolumnData => p_rcolumnData%p_rnextColumn
      call ctab_swapColumns(rtable,p_rcolumnData,p_rcolumnConvRate)      
    else
      ! Throw error if column does not exist
      call output_line("Column does not exist!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_evalConvRateByNumber")
      call sys_halt()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_evalConvRateByObject(rcolumnData,rcolumnConvRate,cevalType)

!<description>
    ! This subroutine computes the convergence rate for column
    ! rcolumnData and stores the result in column rcolumnConvRate. The
    ! evaluation type is specified by cevalType.
!</description>

!<input>
    ! Data column
    type(t_column), intent(in) :: rcolumnData

    ! Type of convergence rate evaluation. One of the CTAB_xxx constants.
    ! CTAB_REDUCTION_RATE:      Reduction rate
    ! CTAB_REDUCTION_RATE_LOG2: Logarithm of reduction rate
    integer, intent(in) :: cevalType
!</input>

!<inputoutput>
    ! Column where convergence rates are stored
    type(t_column), intent(inout) :: rcolumnConvRate
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_cell), pointer :: p_rcellData1,p_rcellData2,p_rcell
    real(DP) :: dvalue1,dvalue2
    real(DP) :: fvalue
    integer(I8)  :: i8value
    integer(I16) :: i16value
    integer(I32) :: i32value
    integer(I64) :: i64value

    ! What type of evaluation should be performed
    select case(cevalType)
    case (CTAB_REDUCTION_RATE)

      ! Get first two cells
      p_rcellData1 => rcolumnData%p_rfirstCell
      p_rcellData2 => p_rcellData1
      if (associated(p_rcellData1)) p_rcellData1 => p_rcellData1%p_rnextCell

      ! Iterate over all cells
      do while(associated(p_rcellData1))
        
        ! Get content from first cell
        select case(p_rcellData1%ctype)
        case(ST_SINGLE)
          dvalue1 = real(transfer(p_rcellData1%svalue,fvalue),DP)
        case (ST_DOUBLE)
          dvalue1 = transfer(p_rcellData1%svalue,dvalue1)
        case (ST_INT8)
          dvalue1 = real(transfer(p_rcellData1%svalue,i8value),DP)
        case (ST_INT16)
          dvalue1 = real(transfer(p_rcellData1%svalue,i16value),DP)
        case (ST_INT32)
          dvalue1 = real(transfer(p_rcellData1%svalue,i32value),DP)
        case (ST_INT64)
          dvalue1 = real(transfer(p_rcellData1%svalue,i64value),DP)
        case default
          ! Throw error if data type is not supported
          call output_line("Unsupported data type!",&
              OU_CLASS_ERROR,OU_MODE_STD,"ctab_evalConvRateByObject")
          call sys_halt()
        end select

        ! Get content from second cell
        select case(p_rcellData2%ctype)
        case(ST_SINGLE)
          dvalue2 = real(transfer(p_rcellData2%svalue,fvalue),DP)
        case (ST_DOUBLE)
          dvalue2 = transfer(p_rcellData2%svalue,dvalue2)
        case (ST_INT8)
          dvalue2 = real(transfer(p_rcellData2%svalue,i8value),DP)
        case (ST_INT16)
          dvalue2 = real(transfer(p_rcellData2%svalue,i16value),DP)
        case (ST_INT32)
          dvalue2 = real(transfer(p_rcellData2%svalue,i32value),DP)
        case (ST_INT64)
          dvalue2 = real(transfer(p_rcellData2%svalue,i64value),DP)
        case default
          ! Throw error if data type is not supported
          call output_line("Unsupported data type!",&
              OU_CLASS_ERROR,OU_MODE_STD,"ctab_evalConvRateByObject")
          call sys_halt()
        end select
        
        ! Add new cell and transfer value
        allocate(p_rcell)
        p_rcell%svalue = transfer(dvalue2/dvalue1,p_rcell%svalue)
        p_rcell%ctype = ST_DOUBLE
        call ctab_addCell(rcolumnConvRate,p_rcell)
        
        ! Proceed with next pair of cells
        p_rcellData2 => p_rcellData1
        p_rcellData1 => p_rcellData1%p_rnextCell
      end do
      
    case (CTAB_REDUCTION_RATE_LOG2)

      p_rcellData1 => rcolumnData%p_rfirstCell
      p_rcellData2 => p_rcellData1
      if (associated(p_rcellData1)) p_rcellData1 => p_rcellData1%p_rnextCell

      do while(associated(p_rcellData1))

        ! Get content from first cell
        select case(p_rcellData1%ctype)
        case(ST_SINGLE)
          dvalue1 = real(transfer(p_rcellData1%svalue,fvalue),DP)
        case (ST_DOUBLE)
          dvalue1 = transfer(p_rcellData1%svalue,dvalue1)
        case (ST_INT8)
          dvalue1 = real(transfer(p_rcellData1%svalue,i8value),DP)
        case (ST_INT16)
          dvalue1 = real(transfer(p_rcellData1%svalue,i16value),DP)
        case (ST_INT32)
          dvalue1 = real(transfer(p_rcellData1%svalue,i32value),DP)
        case (ST_INT64)
          dvalue1 = real(transfer(p_rcellData1%svalue,i64value),DP)
        case default
          ! Throw error if data type is not supported
          call output_line("Unsupported data type!",&
              OU_CLASS_ERROR,OU_MODE_STD,"ctab_evalConvRateByObject")
          call sys_halt()
        end select

        ! Get content from second cell
        select case(p_rcellData2%ctype)
        case(ST_SINGLE)
          dvalue2 = real(transfer(p_rcellData2%svalue,fvalue),DP)
        case (ST_DOUBLE)
          dvalue2 = transfer(p_rcellData2%svalue,dvalue2)
        case (ST_INT8)
          dvalue2 = real(transfer(p_rcellData2%svalue,i8value),DP)
        case (ST_INT16)
          dvalue2 = real(transfer(p_rcellData2%svalue,i16value),DP)
        case (ST_INT32)
          dvalue2 = real(transfer(p_rcellData2%svalue,i32value),DP)
        case (ST_INT64)
          dvalue2 = real(transfer(p_rcellData2%svalue,i64value),DP)
        case default
          ! Throw error if data type is not supported
          call output_line("Unsupported data type!",&
              OU_CLASS_ERROR,OU_MODE_STD,"ctab_evalConvRateByObject")
          call sys_halt()
        end select
        
        ! Add new cell and transfer value
        allocate(p_rcell)
        p_rcell%svalue = transfer(log(dvalue2/dvalue1),p_rcell%svalue)
        p_rcell%ctype = ST_DOUBLE
        call ctab_addCell(rcolumnConvRate,p_rcell)
        
        ! Proceed with next pair of cells
        p_rcellData2 => p_rcellData1
        p_rcellData1 => p_rcellData1%p_rnextCell
      end do

    case default
      call output_line("Unsupported type of evaluation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"ctab_evalConvRateByObject")
      call sys_halt()
    end select

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_done(rtable)

!<description>
    ! This subroutine destroys a convergence table
!</description>

!<inputoutput>
    ! A convergence table
    type(t_convergenceTable), intent(inout) :: rtable
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_column), pointer :: p_rcolumn,p_rcolumn0
    type(t_cell), pointer :: p_rcell,p_rcell0

    rtable%stexCaption = ''
    rtable%stexLabel = ''

    ! Loop through all column
    p_rcolumn => rtable%p_rfirstColumn
    do while (associated(p_rcolumn))

      ! Loop through all cells of the column
      p_rcell => p_rcolumn%p_rfirstCell
      do while (associated(p_rcell))
        
        ! Get next cell
        p_rcell0 => p_rcell
        p_rcell  => p_rcell%p_rnextCell

        ! Destroy current cell
        deallocate(p_rcell0)
      end do

      ! Get next column
      p_rcolumn0 => p_rcolumn
      p_rcolumn  => p_rcolumn%p_rnextColumn

      ! Destroy current column
      deallocate(p_rcolumn0)
    end do

  end subroutine

  !************************************************************************
  ! PRIVATE AUXILIARY ROUTINES
  !************************************************************************

!<subroutine>

  subroutine ctab_addColumn(rtable,rcolumn)

!<description>
    ! This subroutine adds column rcolumn to the end of the table rtable
!</description>

!<inputoutput>
    ! Convergence table
    type(t_convergenceTable), intent(inout) :: rtable

    ! Column object
    type(t_column), intent(inout), target :: rcolumn
!</inputoutput>

!</subroutine>

    if (associated(rtable%p_rfirstColumn)) then
      rtable%p_rlastColumn%p_rnextColumn => rcolumn
      rcolumn%p_rprevColumn => rtable%p_rlastColumn
      rtable%p_rlastColumn  => rcolumn
      rcolumn%p_rnextColumn => null()
    else
      rtable%p_rfirstColumn => rcolumn
      rtable%p_rlastColumn  => rcolumn
      rcolumn%p_rnextColumn => null()
      rcolumn%p_rprevColumn => null()
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ctab_addCell(rcolumn,rcell)

!<description>
    ! This subroutine adds cell rcell to the end of the column rcolumnt
!</description>

!<inputoutput>
    ! Column of the convergence table
    type(t_column), intent(inout) :: rcolumn

    ! Cell object
    type(t_cell), intent(inout), target :: rcell
!</inputoutput>

!</subroutine>

    if (associated(rcolumn%p_rfirstCell)) then
      rcolumn%p_rlastCell%p_rnextCell => rcell
      rcell%p_rprevCell   => rcolumn%p_rlastCell
      rcolumn%p_rlastCell => rcell
      rcell%p_rnextCell   => null()
    else
      rcolumn%p_rfirstCell => rcell
      rcolumn%p_rlastCell  => rcell
      rcell%p_rnextCell    => null()
      rcell%p_rprevCell    => null()
    end if

  end subroutine

  !************************************************************************

!<function>

  function ctab_getNColumns(rtable) result(iresult)

!<description>
    ! This function computes the number of columns of the table
!</description>

!<input>
    ! Convergence table
    type(t_convergenceTable), intent(in) :: rtable
!</input>

!<result>
    ! Number of columns
    integer :: iresult
!</result>

!</function>

    ! local variables
    type(t_column), pointer :: p_rcolumn

    iresult = 0
    p_rcolumn => rtable%p_rfirstColumn

    do while(associated(p_rcolumn))
      iresult = iresult+1
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do

  end function

  !************************************************************************

!<function>

  function ctab_getNRows(rtable) result(iresult)

!<description>
    ! This function computes the number of rows of the table
!</description>

!<input>
    ! Convergence table
    type(t_convergenceTable), intent(in) :: rtable
!</input>

!<result>
    ! Number of rows
    integer :: iresult
!</result>

!</function>

    ! local variables
    type(t_column), pointer :: p_rcolumn

    iresult = 0
    p_rcolumn => rtable%p_rfirstColumn

    do while(associated(p_rcolumn))
      iresult = max(iresult,ctab_getNCellsOfColumn(p_rcolumn))
      p_rcolumn => p_rcolumn%p_rnextColumn
    end do

  end function

  !************************************************************************

!<function>

  function ctab_getNCellsOfColumn(rcolumn) result(iresult)

!<description>
    ! This function computes the number of cells of the column
!</description>

!<input>
    ! Column
    type(t_column), intent(in) :: rcolumn
!</input>

!<result>
    ! Number of cells
    integer :: iresult
!</result>

!</function>

    ! local variables
    type(t_cell), pointer :: p_rcell

    iresult = 0
    p_rcell => rcolumn%p_rfirstCell

    do while(associated(p_rcell))
      iresult = iresult+1
      p_rcell => p_rcell%p_rnextCell
    end do

  end function

  !************************************************************************

!<function>

  function ctab_getColumnWidth(rcolumn) result(iresult)

!<description>
    ! This function computes the width of the column
!</description>

!<input>
    ! Column
    type(t_column), intent(in) :: rcolumn
!</input>

!<result>
    ! Width of the column
    integer :: iresult
!</result>

!</function>

    ! local variables
    type(t_cell), pointer :: p_rcell
    real(SP) :: fvalue
    real(DP) :: dvalue
    integer(I8)  :: i8value
    integer(I16) :: i16value
    integer(I32) :: i32value
    integer(I64) :: i64value

    iresult = len(trim(adjustl(rcolumn%skey)))
    p_rcell => rcolumn%p_rfirstCell

    do while(associated(p_rcell))

      ! What type of data are we?
      select case(p_rcell%ctype)
        
      case (ST_SINGLE)
        fvalue = transfer(p_rcell%svalue,fvalue)
        if (rcolumn%bscientific) then
          iresult = max(iresult,&
                        len(trim(adjustl(sys_sdE(real(fvalue,DP),rcolumn%idigits)))))
        else
          iresult = max(iresult,&
                        len(trim(adjustl(sys_sd(real(fvalue,DP),rcolumn%idigits)))))
        end if
        
      case (ST_DOUBLE)
        dvalue = transfer(p_rcell%svalue,dvalue)
        if (rcolumn%bscientific) then
          iresult = max(iresult,&
                        len(trim(adjustl(sys_sdE(dvalue,rcolumn%idigits)))))
        else
          iresult = max(iresult,&
                        len(trim(adjustl(sys_sd(dvalue,rcolumn%idigits)))))
        end if

      case (ST_INT8)
        i8value = transfer(p_rcell%svalue,i8value)
        iresult = max(iresult,len(trim(adjustl(sys_si(int(i8value),16)))))

      case (ST_INT16)
        i16value = transfer(p_rcell%svalue,i16value)
        iresult = max(iresult,len(trim(adjustl(sys_si(int(i16value),16)))))

      case (ST_INT32)
        i32value = transfer(p_rcell%svalue,i32value)
        iresult = max(iresult,len(trim(adjustl(sys_si(int(i32value),16)))))

      case (ST_INT64)
        i64value = transfer(p_rcell%svalue,i64value)
        iresult = max(iresult,len(trim(adjustl(sys_sli(i64value,16)))))

      case default
        ! Throw error if data type is not supported
        call output_line("Unsupported data type!",&
            OU_CLASS_ERROR,OU_MODE_STD,"ctab_getColumnWidth")
        call sys_halt()
      end select
      
      ! Proceed with next row
      p_rcell => p_rcell%p_rnextCell
    end do

  end function

end module convergencetable

