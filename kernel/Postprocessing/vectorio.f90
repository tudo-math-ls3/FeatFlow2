!#########################################################################
!# ***********************************************************************
!# <name> vectorio </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines and constant necessary to output
!# vectors/arrays to files or read them from files.
!#
!# The following routines can be found in this module:
!#
!# 1.) vecio_writeArray
!#     -> Writes an array into a (text or binary) file
!#
!# 2.) vecio_readArray
!#     -> Reads an array from a (text or binary) file
!#
!# 3.) vecio_writeBlockVectorHR
!#     -> Writes a block vector into a (text or binary) file
!#
!# 4.) vecio_writeVectorHR
!#     -> Writes a scalar vector into a (text or binary) file
!#
!# 5.) vecio_readBlockVectorHR
!#     -> Reads a block vector from a (text or binary) file
!#
!# 6.) vecio_readVectorHR
!#     -> Reads a scalar vector from a (text or binary) file
!#
!# 7.) vecio_writeVectorMaple
!#     -> Writes a scalar vector into a text file in Maple syntax
!#
!# 8.) vecio_writeBlockVectorMaple
!#     -> Writes a block vector into a text file in Maple syntax
!#
!# 9.) vecio_spyVector
!#     -> Writes a scalar vector into a file which can be
!#        visualised by means of the MATLAB command SPY
!#
!# 10.) vecio_spyBlockVector
!#      -> Writes a scalar vector into a file which can be
!#         visualised by means of the MATLAB command SPY
!# </purpose>
!#########################################################################

module vectorio

!$ use omp_lib
  use fsystem
  use genoutput
  use storage
  use io
  use linearsystemscalar
  use linearsystemblock
  use sortstrategybase

  implicit none

  private

  public :: vecio_writeArray
  public :: vecio_readArray
  public :: vecio_writeBlockVectorHR
  public :: vecio_writeVectorHR
  public :: vecio_readBlockVectorHR
  public :: vecio_readVectorHR
  public :: vecio_writeVectorMaple
  public :: vecio_writeBlockVectorMaple
  public :: vecio_spyVector
  public :: vecio_spyBlockVector
  
  interface vecio_writeArray
    module procedure vecio_writeArray_DP
    module procedure vecio_writeArray_SP
    module procedure vecio_writeArray_int
  end interface

  interface vecio_readArray
    module procedure vecio_readArray_DP
    module procedure vecio_readArray_SP
    module procedure vecio_readArray_int
  end interface

contains

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeArray_DP (Ddata, ifile, sfile, sformat, Ipermutation)

  !<description>
    ! Write double precision vector into a text file.
    ! The array is written out 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>

  !<input>
    ! vector: array [:] of double
    real(DP), dimension(:), intent(in) :: Ddata

    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: Permutation for unsorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    integer, dimension(:), optional :: Ipermutation
  !</input>

!</subroutine>

    !local variables
    integer :: i, cf
    real(DP) :: dval

    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE, bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ("Could not open file "//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,"vecio_writeArray_DP")
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (size(Ddata) .le. 0) return

    ! Write the vector.
    ! Unsort the vector on the fly if necessary.
    if (present(sformat)) then
      if (present(Ipermutation)) then
        do i=1, size(Ddata)
          dval = Ddata(Ipermutation(i))
          write (cf,sformat) dval
        end do
      else
        do i=1, size(Ddata)
          dval = Ddata(i)
          write (cf,sformat) dval
        end do
      end if
    else
      if (present(Ipermutation)) then
        do i=1, size(Ddata)
          dval = Ddata(Ipermutation(i))
          write (cf) dval
        end do
      else
        do i=1, size(Ddata)
          dval = Ddata(i)
          write (cf) dval
        end do
      end if
    end if

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeArray_int (Idata, ifile, sfile, sformat, Ipermutation)

  !<description>
    ! Write integer vector into a text file.
    ! The array is written out 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>

  !<input>
    ! vector: array [:] of integer
    integer, dimension(:), intent(in) :: Idata

    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! OPTIONAL: Format string to use for the output; e.g. '(I10)'.
    ! If not specified, data is written to the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: Permutation for unsorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    integer, dimension(:), optional :: Ipermutation
  !</input>

!</subroutine>

    !local variables
    integer :: i, cf
    integer :: ival

    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE, bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeArray_int')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (size(Idata) .le. 0) return

    ! Write the vector.
    ! Unsort the vector on the fly if necessary.
    if (present(sformat)) then
      if (present(Ipermutation)) then
        do i=1, size(Idata)
          ival = Idata(Ipermutation(i))
          write (cf,sformat) ival
        end do
      else
        do i=1, size(Idata)
          ival = Idata(i)
          write (cf,sformat) ival
        end do
      end if
    else
      if (present(Ipermutation)) then
        do i=1, size(Idata)
          ival = Idata(Ipermutation(i))
          write (cf) ival
        end do
      else
        do i=1, size(Idata)
          ival = Idata(i)
          write (cf) ival
        end do
      end if
    end if

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

    ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeArray_SP (Fdata, ifile, sfile, sformat, Ipermutation)

  !<description>
    ! Write single precision vector into a text file.
    ! The array is written out 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>

  !<input>
    ! vector: array [:] of single
    real(SP), dimension(:), intent(in) :: Fdata

    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: Permutation for unsorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    integer, dimension(:), optional :: Ipermutation
  !</input>

!</subroutine>

    !local variables
    integer :: i, cf
    real(SP) :: fval

    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE, bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ("Could not open file "//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,"vecio_writeArray_SP")
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (size(Fdata) .le. 0) return

    ! Write the vector.
    ! Unsort the vector on the fly if necessary.
    if (present(sformat)) then
      if (present(Ipermutation)) then
        do i=1, size(Fdata)
          fval = Fdata(Ipermutation(i))
          write (cf,sformat) fval
        end do
      else
        do i=1, size(Fdata)
          fval = Fdata(i)
          write (cf,sformat) fval
        end do
      end if
    else
      if (present(Ipermutation)) then
        do i=1, size(Fdata)
          fval = Fdata(Ipermutation(i))
          write (cf) fval
        end do
      else
        do i=1, size(Fdata)
          fval = Fdata(i)
          write (cf) fval
        end do
      end if
    end if

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_readArray_DP (Ddata, ifile, sfile, sformat, Ipermutation)

  !<description>
    ! Reads a double precision vector from a file.
    ! The array is read in 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>

  !<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! OPTIONAL: Format string to use for the input; e.g. '(E20.10)'.
    ! If not specified, data is read from the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    ! When reading an array written out by vecio_writeArray_DP,
    ! the format string shall match the setting of the
    ! format string used there.
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: Permutation for sorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    integer, dimension(:), optional :: Ipermutation
  !</input>

  !<output>
    ! Array where to write the data to.
    real(DP), dimension(:), intent(out) :: Ddata
  !</output>

!</subroutine>

    ! local variables
    integer :: i, cf
    real(DP) :: dval

    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, cf, bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ("Could not open file "//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,"vecio_readArray_DP")
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Read the array.
    if (present(sformat)) then
      ! Unsort the vector on the fly if necessary.
      if (present(Ipermutation)) then
        do i=1, size(Ddata)
          read (cf,sformat) dval
          Ddata(Ipermutation(i)) = dval
        end do
      else
        do i=1, size(Ddata)
          read (cf,sformat) dval
          Ddata(i) = dval
        end do
      end if
    else
      ! Unsort the vector on the fly if necessary.
      if (present(Ipermutation)) then
        do i=1, size(Ddata)
          read (cf) dval
          Ddata(Ipermutation(i)) = dval
        end do
      else
        do i=1, size(Ddata)
          read (cf) dval
          Ddata(i) = dval
        end do
      end if
    end if

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_readArray_int (Idata, ifile, sfile, sformat, Ipermutation)

  !<description>
    ! Reads a integer vector from a file.
    ! The array is read in 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>

  !<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! OPTIONAL: Format string to use for the input; e.g. '(E20.10)'.
    ! If not specified, data is read from the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    ! When reading an array written out by vecio_writeArray_DP,
    ! the format string shall match the setting of the
    ! format string used there.
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: Permutation for sorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    integer, dimension(:), optional :: Ipermutation
  !</input>

  !<output>
    ! Array where to write the data to.
    integer, dimension(:), intent(out) :: Idata
  !</output>

!</subroutine>

    ! local variables
    integer :: i, cf
    integer :: ival

    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, cf, bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ("Could not open file "//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,"vecio_readArray_int")
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Read the array.
    if (present(sformat)) then
      ! Unsort the vector on the fly if necessary.
      if (present(Ipermutation)) then
        do i=1, size(Idata)
          read (cf,sformat) ival
          Idata(Ipermutation(i)) = ival
        end do
      else
        do i=1, size(Idata)
          read (cf,sformat) ival
          Idata(i) = ival
        end do
      end if
    else
      ! Unsort the vector on the fly if necessary.
      if (present(Ipermutation)) then
        do i=1, size(Idata)
          read (cf) ival
          Idata(Ipermutation(i)) = ival
        end do
      else
        do i=1, size(Idata)
          read (cf) ival
          Idata(i) = ival
        end do
      end if
    end if

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_readArray_SP (Fdata, ifile, sfile, sformat, Ipermutation)

  !<description>
    ! Reads a single precision vector from a file.
    ! The array is read in 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>

  !<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! OPTIONAL: Format string to use for the input; e.g. '(E20.10)'.
    ! If not specified, data is read from the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    ! When reading an array written out by vecio_writeArray_DP,
    ! the format string shall match the setting of the
    ! format string used there.
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: Permutation for sorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    integer, dimension(:), optional :: Ipermutation
  !</input>

  !<output>
    ! Array where to write the data to.
    real(SP), dimension(:), intent(out) :: Fdata
  !</output>

!</subroutine>

    ! local variables
    integer :: i, cf
    real(SP) :: fval

    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, cf, bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_readArray_SP')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Read the array.
    if (present(sformat)) then
      ! Unsort the vector on the fly if necessary.
      if (present(Ipermutation)) then
        do i=1, size(Fdata)
          read (cf,sformat) fval
          Fdata(Ipermutation(i)) = fval
        end do
      else
        do i=1, size(Fdata)
          read (cf,sformat) fval
          Fdata(i) = fval
        end do
      end if
    else
      ! Unsort the vector on the fly if necessary.
      if (present(Ipermutation)) then
        do i=1, size(Fdata)
          read (cf) fval
          Fdata(Ipermutation(i)) = fval
        end do
      else
        do i=1, size(Fdata)
          read (cf) fval
          Fdata(i) = fval
        end do
      end if
    end if

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeVectorHR (rvector, sarray, bunsort,&
                                  ifile, sfile, sformat, scomment)

  !<description>
    ! This routine writes a (scalar) vector into a (text or binary) file.
  !</description>

  !<input>
    ! The vector to be written out
    type(t_vectorScalar), intent(in) :: rvector

    ! Name of the vector
    character(len=*), intent(in) :: sarray

    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it is unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(in) :: bunsort

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: An additional comment which is added to the 1st line of the file.
    ! Can only be used for formatted files, i.e. if sformat <> "".
    character(len=*), intent(in), optional :: scomment
  !</input>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Ipermutation
    integer :: cf,nchar

    character(len=128) :: S
    character(len=15) :: sarrayname
    character(len=6) :: sformatChar

    if (rvector%NEQ .eq. 0) return ! nothing to do

    ! Open the file
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE,bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeVectorHR')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (present(sformat)) then
      ! Probably write the comment
      if (present(scomment)) then
        write (cf,'(A)') "## "//scomment
      end if

      ! Get length of output strings
      S(:) = ' '
      write (S,sformat) 0.0_DP
      nchar = len(trim(S))

      ! Build array format string
      sformatChar = '(A'//trim(sys_i3(nchar))//')'

      ! Write all format strings into the file
      write (cf,'(A,3A15,2I15)') '# ',sarray, sformat, sformatChar, &
                                  nchar, rvector%NEQ
    else
      sarrayname = sarray
      write (cf) sarrayname, rvector%NEQ
    end if

    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      ! Permuted?
      nullify(p_Ipermutation)
      if (bunsort .and. rvector%bisSorted) then
        ! Get the information how to get the sorted from the unsorted vector.
        call sstrat_getSortedPosInfo (rvector%p_rsortStrategy,p_Ipermutation)
      end if

      call lsyssc_getbase_double (rvector,p_Ddata)

      if (present(sformat)) then
        if (.not. associated(p_Ipermutation)) then
          call vecio_writeArray_DP (p_Ddata, cf, sfile, sformat)
        else
          call vecio_writeArray_DP (p_Ddata, cf, sfile, sformat, p_Ipermutation)
        end if
      else
        if (.not. associated(p_Ipermutation)) then
          call vecio_writeArray_DP (p_Ddata, cf, sfile)
        else
          call vecio_writeArray_DP (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
        end if
      end if
    case default
      call output_line ('Unsupported vector precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeVectorHR')
      call sys_halt()
    end select

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_readVectorHR (rvector, sarray, bunsort,&
                                 ifile, sfile, bformatted)

  !<description>
    ! This routine reads a (scalar) vector from a (text or binary) file.
    !
    ! Note: The input data in the file may be written out with vecio_writeVectorHR
    ! or vecio_writeBlockVectorHR; in the latter case, a preciously written
    ! out block vector is read in as scalar vector.
  !</description>

  !<input>
    ! Input channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Read from channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! Name of the file where to read from. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! Read unsorted vector.
    ! =TRUE:  Read data and sort it according to the sorting strategy
    !         in the vector (if the vector is sorted)
    ! =FALSE: Read vector as it is.
    logical, intent(in) :: bunsort

    ! Whether to read data formatted or unformatted.
    ! TRUE  = Treat data in input file formatted, i.e. in human readable form.
    ! FALSE = Data in the input file is unformatted, i.e. in processor
    !         dependent form.
    ! A vector written out by vecio_writeVectorHR with a format specifier
    ! sformat being specified shall be read with bformatted=TRUE.
    ! A vector written out by vecio_writeVectorHR without a format specifier
    ! sformat being specified shall be read with bformatted=FALSE.
    logical, intent(in) :: bformatted
  !</input>

  !<inputoutput>
    ! The vector to be read in.
    ! If the vector is not initialised, a new vector is automatically created
    ! with the correct size.
    type(t_vectorScalar), intent(inout) :: rvector
  !</inputoutput>

  !<output>
    ! Name of the vector
    character(len=*), intent(out) :: sarray
  !</output>


!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Ipermutation
    integer :: cf,nchar
    integer :: NEQ

    character(len=15) :: sarrayname,sformat
    character(len=6) :: sformatChar,S

    ! Open the file
    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, cf, bformatted=bformatted)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeVectorHR')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (bformatted) then
      ! Peek the first line(s). Ignore all comments.
      do
        read (cf,'(A2)',ADVANCE="NO") S
        if (S .eq. "##") then
          read (cf,*)
        else
          exit
        end if
      end do

      ! Get the format specification from the file
      read (cf,'(3A15,2I15)') sarrayname, sformat, sformatChar, &
                                nchar, NEQ
    else
      ! Get array information from the file
      read (cf) sarrayname, NEQ
    end if

    sarray = sarrayname

    ! Does the vector exist? If not, we create a new one.
    if (rvector%NEQ .eq. 0) then
      call lsyssc_createVector (rvector,NEQ,.false.,ST_DOUBLE)
    end if

    if (rvector%NEQ .ne. NEQ) then
      call output_line ('Vector has wrong size', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeVectorHR')
      call sys_halt()
    end if

    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      ! Permuted?
      nullify(p_Ipermutation)
      if (bunsort .and. rvector%bisSorted) then
        ! Get the information how to get the sorted from the unsorted vector.
        call sstrat_getSortedPosInfo (rvector%p_rsortStrategy,p_Ipermutation)
      end if

      call lsyssc_getbase_double (rvector,p_Ddata)

      if (bformatted) then
        if (.not. associated(p_Ipermutation)) then
          call vecio_readArray_DP (p_Ddata, cf, sfile, sformat)
        else
          call vecio_readArray_DP (p_Ddata, cf, sfile, sformat, p_Ipermutation)
        end if
      else
        if (.not. associated(p_Ipermutation)) then
          call vecio_readArray_DP (p_Ddata, cf, sfile)
        else
          call vecio_readArray_DP (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
        end if
      end if
    case default
      call output_line ('Unsupported vector precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_readVectorHR')
      call sys_halt()
    end select

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine vecio_readVectorHR

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeBlockVectorHR (rvector, sarray, bunsort,&
                                       ifile, sfile, sformat,scomment)

  !<description>
    ! This routine writes a block vector into a (text or binary) file.
    ! The output file can be read in with vecio_readBlockVectorHR as
    ! block vector or with vecio_readVectorHR as scalar vector.
  !</description>

  !<input>
    ! The vector to be written out
    type(t_vectorBlock), intent(in) :: rvector

    ! Name of the vector
    character(len=*), intent(in) :: sarray

    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it is unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(in) :: bunsort

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(in), optional :: sformat

    ! OPTIONAL: An additional comment which is added to the 1st line of the file.
    ! Can only be used for formatted files, i.e. if sformat <> "".
    character(len=*), intent(in), optional :: scomment
  !</input>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Ipermutation
    integer :: cf,nchar,I

    character(len=128) :: S
    character(len=15) :: sarrayname
    character(len=6) :: sformatChar

    if (rvector%NEQ .eq. 0) return ! nothing to do

    ! Open the file
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE,bformatted=present(sformat))
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeBlockVectorHR')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Write all format strings into the file
    if (present(sformat)) then
      ! Probably write the comment
      if (present(scomment)) then
        write (cf,'(A)') "## "//scomment
      end if

      ! Get length of output strings
      S(:) = ' '
      write (S,sformat) 0.0_DP
      nchar = len(trim(S))

      ! Build array format string
      sformatChar = '(A'//trim(sys_i3(nchar))//')'

      write (cf,'(A,3A15,3I15)',ADVANCE="NO") '# ',sarray, sformat, sformatChar, &
                                               nchar, rvector%NEQ, rvector%nblocks
      ! Write block structure
      do i=1,rvector%nblocks
        write (cf,'(I15)',ADVANCE="NO") rvector%RvectorBlock(i)%NEQ
      end do

      ! New line
      write (cf,*)

    else
      sarrayname = sarray
      write (cf) sarrayname, rvector%NEQ,rvector%nblocks
      ! Write block structure
      write (cf) (rvector%RvectorBlock(i)%NEQ,i=1,rvector%nblocks)
    end if

    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      do i=1,rvector%nblocks
        ! Permuted?
        nullify(p_Ipermutation)
        if (bunsort .and. rvector%RvectorBlock(i)%bisSorted) then
          ! Get the information how to get the sorted from the unsorted vector.
          call sstrat_getSortedPosInfo (rvector%RvectorBlock(i)%p_rsortStrategy,p_Ipermutation)
        end if

        call lsyssc_getbase_double (rvector%RvectorBlock(i),p_Ddata)

        if (present(sformat)) then
          if (.not. associated(p_Ipermutation)) then
            call vecio_writeArray_DP (p_Ddata, cf, sfile, sformat)
          else
            call vecio_writeArray_DP (p_Ddata, cf, sfile, sformat, p_Ipermutation)
          end if
        else
          if (.not. associated(p_Ipermutation)) then
            call vecio_writeArray_DP (p_Ddata, cf, sfile)
          else
            call vecio_writeArray_DP (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
          end if
        end if
      end do
    case default
      call output_line ('Unsupported vector precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeBlockVectorHR')
      call sys_halt()
    end select

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_readBlockVectorHR (rvector, sarray, bunsorted,&
                                      ifile, sfile, bformatted)

  !<description>
    ! This routine reads a block vector from a (text or binary) file.
    !
    ! rvector may not be initialised when calling this routine; in this
    ! case, if the source file specifies a block vector, a new block
    ! vector is created in the same structure as written out
    ! by vecio_writeBlockVectorHR.
    !
    ! If rvector is initialised, the data is read from the file according
    ! to the structure of rvector.
    !
    ! If rvector specifies a structure and bformatted=TRUE, this routine
    ! is compatible to the FEAT 1.0 format of formatted vectors!
  !</description>

  !<input>
    ! Input channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Read from channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! Name of the file where to read from. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! Read unsorted vector.
    ! =TRUE:  Read data and sort it according to the sorting strategy
    !         in the vector (if the vector is sorted)
    ! =FALSE: Read vector as it is.
    logical, intent(in) :: bunsorted

    ! Whether to read data formatted or unformatted.
    ! TRUE  = Treat data in input file formatted, i.e. in human readable form.
    ! FALSE = Data in the input file is unformatted, i.e. in processor
    !         dependent form.
    ! A vector written out by vecio_writeBlockVectorHR with a format specifier
    ! sformat being specified shall be read with bformatted=TRUE.
    ! A vector written out by vecio_writeBlockVectorHR without a format specifier
    ! sformat being specified shall be read with bformatted=FALSE.
    logical, intent(in) :: bformatted
  !</input>

  !<inputoutput>
    ! The vector to be read in.
    ! If the vector is not initialised, a new vector is automatically created
    ! with the correct size.
    type(t_vectorBlock), intent(inout) :: rvector
  !</inputoutput>

  !<output>
    ! Name of the vector
    character(len=*), intent(out) :: sarray
  !</output>


!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Ipermutation
    integer :: cf,nchar,nblocks,i
    integer :: NEQ
    integer, dimension(:), allocatable :: IblockSize

    character(len=15) :: sarrayname,sformat
    character(len=6) :: sformatChar,S

    ! Open the file
    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, cf, bformatted=bformatted)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_readBlockVectorHR')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (bformatted) then
      ! Peek the first line(s). Ignore all comments.
      do
        read (cf,'(A2)',ADVANCE="NO") S
        if (S .eq. "##") then
          read (cf,*)
        else
          exit
        end if
      end do

      ! Get the format specification from the file
      read (cf,'(3A15,3I15)',ADVANCE="NO") sarrayname, sformat, sformatChar, &
                                             nchar, NEQ
      ! If rvector does not specify the structure, read the vector structure
      ! from the file.
      if (rvector%NEQ .eq. 0) then
        read (cf,'(I15)',ADVANCE="NO") nblocks

        ! Get the block structure
        allocate(IblockSize(nblocks))
        do i=1,nblocks
          read (cf,'(I15)',ADVANCE="NO") IblockSize(i)
        end do

        ! New line
        read (cf,*)
      else
        ! Take the structure from rvector
        nblocks = rvector%nblocks
        allocate(IblockSize(nblocks))
        IblockSize(:) = rvector%RvectorBlock(1:nblocks)%NEQ

        ! New line
        read (cf,*)
      end if
    else
      ! Get array information from the file
      read (cf) sarrayname, NEQ, nblocks

      ! Get the block structure
      allocate(IblockSize(nblocks))
      read (cf) (IblockSize(i),i=1,nblocks)
    end if

    sarray = sarrayname

    ! Does the vector exist? If not, we create a new one.
    if (rvector%NEQ .eq. 0) then
      call lsysbl_createVecBlockDirect (rvector,IblockSize,.false.,ST_DOUBLE)
    end if

    ! Size of vector must match! Size of subvectors not -- it is not a bug,
    ! it is a feature ;-)
    if (rvector%NEQ .ne. NEQ) then
      call output_line ('Vector has wrong size!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_readBlockVectorHR')
      call sys_halt()
    end if

    ! We do not need the block size anymore.
    deallocate (IblockSize)

    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      do i=1,nblocks
        ! Permuted?
        nullify(p_Ipermutation)
        if (bunsorted .and. rvector%RvectorBlock(i)%bisSorted) then
          ! Get the information how to get the sorted from the unsorted vector.
          call sstrat_getSortedPosInfo (rvector%RvectorBlock(i)%p_rsortStrategy,p_Ipermutation)
        end if

        call lsyssc_getbase_double (rvector%RvectorBlock(i),p_Ddata)

        if (bformatted) then
          if (.not. associated(p_Ipermutation)) then
            call vecio_readArray_DP (p_Ddata, cf, sfile, sformat)
          else
            call vecio_readArray_DP (p_Ddata, cf, sfile, sformat, p_Ipermutation)
          end if
        else
          if (.not. associated(p_Ipermutation)) then
            call vecio_readArray_DP (p_Ddata, cf, sfile)
          else
            call vecio_readArray_DP (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
          end if
        end if
      end do
    case default
      call output_line ('Unsupported vector precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_readBlockVectorHR')
      call sys_halt()
    end select

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeVectorMaple (rvector, sarray, bunsort,&
                                     ifile, sfile, sformat)

  !<description>
    ! This routine writes a (scalar) vector into a file in the MAPLE syntax.
  !</description>

  !<input>
    ! The vector to be written out
    type(t_vectorScalar), intent(in) :: rvector

    ! Name of the vector
    character(len=*), intent(in) :: sarray

    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it is unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(in) :: bunsort

    ! Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(in) :: sformat
  !</input>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Ipermutation
    integer :: cf,nchar

    character(len=128) :: S
    character(len=6) :: sformatChar

    if (rvector%NEQ .eq. 0) return ! nothing to do

    ! Open the file
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE,bformatted=.true.)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeVectorMaple')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Get length of output strings
    S(:) = ' '
    write (S,sformat) 0.0_DP
    nchar = len(trim(S))

    ! Build array format string
    sformatChar = '(A'//trim(sys_i3(nchar))//')'

    ! Write a header:
    write (cf,'(A)',ADVANCE='NO') sarray//':=vector([';

    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      ! Permuted?
      nullify(p_Ipermutation)
      if (bunsort .and. rvector%bisSorted) then
        ! Get the information how to get the sorted from the unsorted vector.
        call sstrat_getSortedPosInfo (rvector%p_rsortStrategy,p_Ipermutation)
      end if

      call lsyssc_getbase_double (rvector,p_Ddata)

      if (.not. associated(p_Ipermutation)) then
        call vecio_writeMapleArray_DP (p_Ddata, cf, sformat)
      else
        call vecio_writeMapleArray_DP (p_Ddata, cf, sformat, p_Ipermutation)
      end if

      ! Footer
      write (cf,'(A)') '):';

    case default
      call output_line ('Unsupported vector precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeVectorMaple')
      call sys_halt()
    end select

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeBlockVectorMaple (rvector, sarray, bunsort,&
                                          ifile, sfile, sformat)

  !<description>
    ! This routine writes a (block) vector into a file in the MAPLE syntax.
  !</description>

  !<input>
    ! The vector to be written out
    type(t_vectorBlock), intent(in) :: rvector

    ! Name of the vector
    character(len=*), intent(in) :: sarray

    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile

    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile

    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it is unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(in) :: bunsort

    ! Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(in) :: sformat
  !</input>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Ipermutation
    integer :: cf,nchar,iblock

    character(len=128) :: S
    character(len=6) :: sformatChar

    if (rvector%NEQ .eq. 0) return ! nothing to do

    ! Open the file
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE,bformatted=.true.)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeBlockVectorMaple')
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Get length of output strings
    S(:) = ' '
    write (S,sformat) 0.0_DP
    nchar = len(trim(S))

    ! Build array format string
    sformatChar = '(A'//trim(sys_i3(nchar))//')'

    ! Write a header:
    write (cf,'(A)',ADVANCE='NO') sarray//':=vector([';

    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      ! Loop over the blocks
      do iblock = 1,rvector%nblocks
        ! Permuted?
        nullify(p_Ipermutation)
        if (bunsort .and. rvector%RvectorBlock(iblock)%bisSorted) then
          ! Get the information how to get the sorted from the unsorted vector.
          call sstrat_getSortedPosInfo (rvector%RvectorBlock(iblock)%p_rsortStrategy,p_Ipermutation)
        end if

        call lsyssc_getbase_double (rvector%RvectorBlock(iblock),p_Ddata)

        if (.not. associated(p_Ipermutation)) then
          call vecio_writeMapleArray_DP (p_Ddata, cf, sformat)
        else
          call vecio_writeMapleArray_DP (p_Ddata, cf, sformat, p_Ipermutation)
        end if

        ! If this is not the last block, attach more data
        if (iblock .ne. rvector%nblocks) then
          write (cf,'(A)', ADVANCE='NO') ','
        end if

      end do
      ! Footer
      write (cf,'(A)') '):';

    case default
      call output_line ('Unsupported vector precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vecio_writeBlockVectorMaple')
      call sys_halt()
    end select

    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecio_writeMapleArray_DP (Ddata, ifile, sformat, Ipermutation)

!<description>
  ! INTERNAL SUBROUTINE.
  ! Writes the data of an array to the Maple output file iodentified by the
  ! output channel ifile.
!</description>

!<input>
  ! vector: array [:] of double
  real(DP), dimension(:), intent(in) :: Ddata

  ! output channel to use for output
  integer, intent(in) :: ifile

  ! Format string to use for the output; e.g. '(E20.10)'.
  ! If not specified, data is written to the file unformatted
  ! (i.e. in a computer dependent, not human readable form).
  character(len=*), intent(in) :: sformat

  ! OPTIONAL: Permutation for unsorting.
  ! If specified, this permutation tells how to unsort a vector before
  ! writing it to the file.
  integer, dimension(:), optional :: Ipermutation
!</input>

!</subroutine>

    !local variables
    integer :: i, cf
    character(LEN=32) :: sdata
    real(DP) :: dval

    cf = ifile

    if (size(Ddata) .le. 0) return

    ! Write the vector.
    ! Unsort the vector on the fly if necessary.
    if (present(Ipermutation)) then
      do i=1, size(Ddata)-1
        dval = Ddata(Ipermutation(i))
        write (sdata,sformat) dval
        write (cf,'(A)',ADVANCE='NO') trim(adjustl(sdata))//','
      end do
      dval = Ddata(Ipermutation(i))
      write (sdata,sformat) dval
      write (cf,'(A)',ADVANCE='NO') trim(adjustl(sdata));
    else
      do i=1, size(Ddata)-1
        dval = Ddata(i)
        write (sdata,sformat) dval
        write (cf,'(A)',ADVANCE='NO') trim(adjustl(sdata))//','
      end do
      dval = Ddata(i)
      write (sdata,sformat) dval
      write (cf,'(A)',ADVANCE='NO') trim(adjustl(sdata));
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecio_spyVector(sfilename, svectorName, rvector,&
      bunsort, cstatus, dthreshold)

!<description>

    ! Writes a scalar vector to a file so that it can be visualised in
    ! MATLAB by means of the SPY command.
    !
    ! To load a vector written out by this routine, one has simply to
    ! type the name of the ".m"-file that is written out. MATLAB will
    ! read the file and create a sparse vector with the name
    ! svectorName in memory.
!</description>

!<input>
    ! File name of the MATLAB file with fileextension.
    character(LEN=*), intent(in) :: sfileName

    ! Name of the vector in MATLAB file. This will be the name of the
    ! variable containing the vector data when reading the file into matlab.
    character(LEN=*), intent(in) :: svectorName

    ! Source vector
    type(t_vectorScalar), intent(in) :: rvector

    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it is unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(in) :: bunsort

    ! OPTIONAL: status of file
    integer, intent(in), optional :: cstatus

    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(in), optional :: dthreshold
!</input>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_Fa
     integer, dimension(:), pointer :: p_Ipermutation
    integer :: iunit,ieq,ivar
    real(DP) :: dthres
    character(LEN=10) :: cstat,cpos

    if (rvector%NEQ .eq. 0) return ! nothing to do

    ! Replace small values by zero
    dthres = 1E-12_DP
    if (present(dthreshold)) dthres = dthreshold

    ! Set file status if required
    if (present(cstatus)) then
      select case(cstatus)
      case (IO_NEW)
        cstat="NEW"; cpos="ASIS"
      case (IO_REPLACE)
        cstat="REPLACE"; cpos="ASIS"
      case (IO_OLD)
        cstat="OLD"; cpos="APPEND"
      case default
        cstat="UNKNOWN"; cpos ="ASIS"
      end select
    else
      cstat="UNKNOWN"; cpos="ASIS"
    end if

    ! Open output file
    iunit=sys_getFreeUnit()
    open (UNIT=iunit,STATUS=trim(cstat),POSITION=trim(cpos),FILE=trim(adjustl(sfilename)))

    ! Which vector type are we?
    select case(rvector%cdataType)
    case (ST_DOUBLE)
      write(UNIT=iunit,FMT=10)

      ! Permuted?
      nullify(p_Ipermutation)
      if (bunsort .and. rvector%bisSorted) then
        ! Get the information how to get the sorted from the unsorted vector.
        call sstrat_getSortedPosInfo (rvector%p_rsortStrategy,p_Ipermutation)
      end if

      call lsyssc_getbase_double(rvector,p_Da)

      if (.not. associated(p_Ipermutation)) then

        do ieq=1,rvector%NEQ
          do ivar=1,rvector%NVAR
            if (abs(p_Da(rvector%NVAR*(ieq-1)+ivar)) .ge. dthres) then
              write(UNIT=iunit,FMT=20) 1,rvector%NVAR*(ieq-1)+ivar,&
                  p_Da(rvector%NVAR*(ieq-1)+ivar)
            end if
          end do
        end do

      else

        do ieq=1,rvector%NEQ
          do ivar=1,rvector%NVAR
            if (abs(p_Da(rvector%NVAR*(ieq-1)+ivar)) .ge. dthres) then
              write(UNIT=iunit,FMT=20) 1,rvector%NVAR*(p_Ipermutation(ieq)-1)+ivar,&
                  p_Da(rvector%NVAR*(p_Ipermutation(ieq)-1)+ivar)
            end if
          end do
        end do

      end if

      write(UNIT=iunit,FMT=30)

    case (ST_SINGLE)
      write(UNIT=iunit,FMT=10)

      ! Permuted?
      nullify(p_Ipermutation)
      if (bunsort .and. rvector%bisSorted) then
        ! Get the information how to get the sorted from the unsorted vector.
        call sstrat_getSortedPosInfo (rvector%p_rsortStrategy,p_Ipermutation)
      end if

      call lsyssc_getbase_single(rvector,p_Fa)

      if (.not. associated(p_Ipermutation)) then

        do ieq=1,rvector%NEQ
          do ivar=1,rvector%NVAR
            if (abs(p_Fa(rvector%NVAR*(ieq-1)+ivar)) .ge. real(dthres,SP)) then
              write(UNIT=iunit,FMT=20) 1,rvector%NVAR*(ieq-1)+ivar,&
                  p_Fa(rvector%NVAR*(ieq-1)+ivar)
            end if
          end do
        end do

      else

        do ieq=1,rvector%NEQ
          do ivar=1,rvector%NVAR
            if (abs(p_Fa(rvector%NVAR*(ieq-1)+ivar)) .ge. real(dthres,SP)) then
              write(UNIT=iunit,FMT=20) 1,rvector%NVAR*(p_Ipermutation(ieq)-1)+ivar,&
                  p_Fa(rvector%NVAR*(p_Ipermutation(ieq)-1)+ivar)
            end if
          end do
        end do

      end if

      write(UNIT=iunit,FMT=30)

    case default
      call output_line ('Unsupported vector precision!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'vecio_spyVector')
      call sys_halt()
    end select

    ! Close file
    write(UNIT=iunit,FMT=40) svectorName, 1, rvector%NEQ*rvector%NVAR,&
                             svectorName, rvector%NEQ*rvector%NVAR, 1
    close(UNIT=iunit)

10  format("data=[...")
20  format(I10,1X,I10,1X,E15.8,";...")
30  format("];")
40  format("if ~isempty(data), ",&
        A,"=transpose(sparse(data(:,1),data(:,2),data(:,3),",I10,",",I10,"));",&
        "else ",&
        A,"=sparse(",I10,",",I10,");",&
        "end; clear data;")
    
  end subroutine vecio_spyVector

  ! ***************************************************************************

!<subroutine>

  subroutine vecio_spyBlockVector(sfilename, svectorName, rvector,&
      bunsort, cstatus, dthreshold)

!<description>

    ! Writes a block vector to a file so that it can be visualised in
    ! MATLAB by means of the SPY command.
    !
    ! To load a vector written out by this routine, one has simply to
    ! type the name of the ".m"-file that is written out. MATLAB will
    ! read the file and create a sparse vector with the name
    ! svectorName in memory.
!</description>

!<input>
    ! File name of the MATLAB file with fileextension.
    character(LEN=*), intent(in) :: sfileName

    ! Name of the vector in MATLAB file. This will be the name of the
    ! variable containing the vector data when reading the file into matlab.
    character(LEN=*), intent(in) :: svectorName

    ! Source vector
    type(t_vectorBlock), intent(in) :: rvector

    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it is unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(in) :: bunsort

    ! OPTIONAL: status of file
    integer, intent(in), optional :: cstatus

    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(in), optional :: dthreshold
!</input>
!</subroutine>

    ! local variables
    integer :: iunit,i
    
    ! Spy scalar subvectors
    do i = 1, rvector%nblocks
      if (i .eq. 1) then
        call vecio_spyVector(sfilename,&
            svectorName//"_"//trim(adjustl(sys_si(i,8))),&
            rvector%RvectorBlock(i), bunsort, cstatus, dthreshold)
      else
        call vecio_spyVector(sfilename,&
            svectorName//"_"//trim(adjustl(sys_si(i,8))),&
            rvector%RvectorBlock(i), bunsort, IO_OLD, dthreshold)
      end if
    end do

    ! Open output file
    iunit=sys_getFreeUnit()
    open (UNIT=iunit,STATUS="OLD",POSITION="APPEND",FILE=trim(adjustl(sfilename)))

    write(UNIT=iunit,FMT=10) svectorName
    do i = 1,rvector%nblocks
      write(UNIT=iunit,FMT=20) svectorName//"_"//trim(adjustl(sys_si(i,8)))
    end do
    write(UNIT=iunit,FMT=30)
    
    ! Close file
    close(UNIT=iunit)
    
10  format(A,"=[...")
20  format(A,";...")
30  format("];")

  end subroutine vecio_spyBlockVector

end module
