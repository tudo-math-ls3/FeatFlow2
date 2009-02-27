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
!# 1.) vecio_writeArray_Dble
!#     -> Writes an array into a (text or binary) file
!#
!# 2.) vecio_readArray_Dble
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
!# 7.) vecio_writeBlockVectorMaple
!#     -> Writes a block vector into a text file in Maple syntax
!# </purpose>
!#########################################################################

module vectorio

  use fsystem
  use storage
  use io
  use linearsystemscalar
  use linearsystemblock
  
  implicit none

  contains

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeArray_Dble (Ddata, ifile, sfile, sformat, Ipermutation)
  
  !<description>
    ! Write double precision vector into a text file.
    ! The array is written out 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>
    
  !<input>
    ! vector: array [:] of double
    real(DP), dimension(:), intent(IN) :: Ddata
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(IN), optional :: sformat

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
        print *, 'vecio_writeArray_Dble: Could not open file '// &
                 trim(sfile)
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
  subroutine vecio_readArray_Dble (Ddata, ifile, sfile, sformat, Ipermutation)
  
  !<description>
    ! Reads a double precision vector from a file.
    ! The array is read in 'as it is', i.e. as specified by sformat
    ! without any additional header or footer.
  !</description>
    
  !<input>
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! OPTIONAL: Format string to use for the input; e.g. '(E20.10)'.
    ! If not specified, data is read from the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    ! When reading an array written out by vecio_writeArray_Dble,
    ! the format string shall match the setting of the 
    ! format string used there.
    character(len=*), intent(IN), optional :: sformat
    
    ! OPTIONAL: Permutation for sorting.
    ! If specified, this permutation tells how to unsort a vector before
    ! writing it to the file.
    integer, dimension(:), optional :: Ipermutation
  !</input>
  
  !<output>
    ! Array where to write the data to.
    real(DP), dimension(:), intent(OUT) :: Ddata
  !</output>
    
!</subroutine>
    
    ! local variables
    integer :: i, cf
    real(DP) :: dval
    
    if (ifile .eq. 0) then
      call io_openFileForReading(sfile, cf, bformatted=present(sformat))
      if (cf .eq. -1) then
        print *, 'vecio_readArray_Dble: Could not open file '// &
                 trim(sfile)
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
  subroutine vecio_writeVectorHR (rvector, sarray, bunsort,&
                                  ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a (scalar) vector into a (text or binary) file.
  !</description>
    
  !<input>
    ! The vector to be written out
    type(t_vectorScalar), intent(IN) :: rvector
    
    ! Name of the vector
    character(len=*), intent(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it's unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(IN) :: bunsort

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(IN), optional :: sformat
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
        print *, 'vecio_writeVectorHR: Could not open file '// &
                 trim(sfile)
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (present(sformat)) then
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
      if (bunsort .and. (lsyssc_isVectorSorted (rvector))) then
        call storage_getbase_int (rvector%h_IsortPermutation,p_Ipermutation)
        ! We must use the inverse permutation
        p_Ipermutation => p_Ipermutation(rvector%NEQ+1:)
      end if

      call lsyssc_getbase_double (rvector,p_Ddata)
      
      if (present(sformat)) then
        if (.not. associated(p_Ipermutation)) then
          call vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat)
        else 
          call vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
        end if
      else
        if (.not. associated(p_Ipermutation)) then
          call vecio_writeArray_Dble (p_Ddata, cf, sfile)
        else 
          call vecio_writeArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
        end if
      end if
    case DEFAULT
      print *,'vecio_writeVectorHR: Unsupported vector precision.'
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
    ! <> 0: Read from channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to read from. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Read unsorted vector.
    ! =TRUE:  Read data and sort it according to the sorting strategy
    !         in the vector (if the vector is sorted)
    ! =FALSE: Read vector as it is.
    logical, intent(IN) :: bunsort

    ! Whether to read data formatted or unformatted.
    ! TRUE  = Treat data in input file formatted, i.e. in human readable form.
    ! FALSE = Data in the input file is unformatted, i.e. in processor
    !         dependent form.
    ! A vector written out by vecio_writeVectorHR with a format specifier
    ! sformat being specified shall be read with bformatted=TRUE.
    ! A vector written out by vecio_writeVectorHR without a format specifier
    ! sformat being specified shall be read with bformatted=FALSE.
    logical, intent(IN) :: bformatted
  !</input>

  !<inputoutput>
    ! The vector to be read in.
    ! If the vector is not initialised, a new vector is automatically created
    ! with the correct size.
    type(t_vectorScalar), intent(INOUT) :: rvector
  !</inputoutput>

  !<output>    
    ! Name of the vector
    character(len=*), intent(OUT) :: sarray
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
        print *, 'vecio_writeVectorHR: Could not open file '// &
                 trim(sfile)
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (bformatted) then
      ! Get the format specification from the file
      read (cf,'(A2,3A15,2I15)') S,sarrayname, sformat, sformatChar, &
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
      print *,'vecio_readVectorHR: Vector has wrong size!'
      call sys_halt()
    end if
    
    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      ! Permuted?
      nullify(p_Ipermutation)
      if (bunsort .and. (lsyssc_isVectorSorted (rvector))) then
        call storage_getbase_int (rvector%h_IsortPermutation,p_Ipermutation)
        ! We must use the inverse permutation
        p_Ipermutation => p_Ipermutation(NEQ+1:)
      end if

      call lsyssc_getbase_double (rvector,p_Ddata)
      
      if (bformatted) then
        if (.not. associated(p_Ipermutation)) then
          call vecio_readArray_Dble (p_Ddata, cf, sfile, sformat)
        else 
          call vecio_readArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
        end if
      else
        if (.not. associated(p_Ipermutation)) then
          call vecio_readArray_Dble (p_Ddata, cf, sfile)
        else 
          call vecio_readArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
        end if
      end if
    case DEFAULT
      print *,'vecio_readVectorHR: Unsupported vector precision.'
      call sys_halt()
    end select
    
    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)
    
  end subroutine 

  ! ***************************************************************************

!<subroutine>
  subroutine vecio_writeBlockVectorHR (rvector, sarray, bunsort,&
                                       ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a block vector into a (text or binary) file.
    ! The output file can be read in with vecio_readBlockVectorHR as
    ! block vector or with vecio_readVectorHR as scalar vector.
  !</description>
    
  !<input>
    ! The vector to be written out
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! Name of the vector
    character(len=*), intent(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it's unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(IN) :: bunsort

    ! OPTIONAL: Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(IN), optional :: sformat
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
        print *, 'vecio_writeBlockVectorHR: Could not open file '// &
                 trim(sfile)
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    ! Write all format strings into the file
    if (present(sformat)) then
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
        if (bunsort .and. (lsyssc_isVectorSorted (rvector%RvectorBlock(i)))) then
          call storage_getbase_int (&
              rvector%RvectorBlock(i)%h_IsortPermutation,p_Ipermutation)
          ! We must use the inverse permutation
          p_Ipermutation => p_Ipermutation(rvector%RvectorBlock(i)%NEQ+1:)
        end if

        call lsyssc_getbase_double (rvector%RvectorBlock(i),p_Ddata)
        
        if (present(sformat)) then
          if (.not. associated(p_Ipermutation)) then
            call vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat)
          else 
            call vecio_writeArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
          end if
        else
          if (.not. associated(p_Ipermutation)) then
            call vecio_writeArray_Dble (p_Ddata, cf, sfile)
          else 
            call vecio_writeArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
          end if
        end if
      end do
    case DEFAULT
      print *,'vecio_writeBlockVectorHR: Unsupported vector precision.'
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
    ! <> 0: Read from channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to read from. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Read unsorted vector.
    ! =TRUE:  Read data and sort it according to the sorting strategy
    !         in the vector (if the vector is sorted)
    ! =FALSE: Read vector as it is.
    logical, intent(IN) :: bunsorted

    ! Whether to read data formatted or unformatted.
    ! TRUE  = Treat data in input file formatted, i.e. in human readable form.
    ! FALSE = Data in the input file is unformatted, i.e. in processor
    !         dependent form.
    ! A vector written out by vecio_writeBlockVectorHR with a format specifier
    ! sformat being specified shall be read with bformatted=TRUE.
    ! A vector written out by vecio_writeBlockVectorHR without a format specifier
    ! sformat being specified shall be read with bformatted=FALSE.
    logical, intent(IN) :: bformatted
  !</input>

  !<inputoutput>
    ! The vector to be read in.
    ! If the vector is not initialised, a new vector is automatically created
    ! with the correct size.
    type(t_vectorBlock), intent(INOUT) :: rvector
  !</inputoutput>

  !<output>    
    ! Name of the vector
    character(len=*), intent(OUT) :: sarray
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
        print *, 'vecio_writeVectorHR: Could not open file '// &
                 trim(sfile)
        call sys_halt()
      end if
    else
      cf = ifile
    end if

    if (bformatted) then
      ! Get the format specification from the file
      read (cf,'(A2,3A15,3I15)',ADVANCE="NO") S,sarrayname, sformat, sformatChar, &
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
    
    ! Size of vector must match! Size of subvectors not -- it's not a bug,
    ! it's a feature ;-)
    if (rvector%NEQ .ne. NEQ) then
      print *,'vecio_readBlockVectorHR: Vector has wrong size!'
      call sys_halt()
    end if

    ! We don't need the block size anymore.    
    deallocate (IblockSize)
    
    ! Vector precision?
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      do i=1,nblocks
        ! Permuted?
        nullify(p_Ipermutation)
        if (bunsorted .and. (lsyssc_isVectorSorted (rvector%RvectorBlock(i)))) then
          call storage_getbase_int (rvector%RvectorBlock(i)%h_IsortPermutation,&
              p_Ipermutation)
          ! We must use the inverse permutation
          p_Ipermutation => p_Ipermutation(rvector%RvectorBlock(i)%NEQ+1:)
        end if

        call lsyssc_getbase_double (rvector%RvectorBlock(i),p_Ddata)
        
        if (bformatted) then
          if (.not. associated(p_Ipermutation)) then
            call vecio_readArray_Dble (p_Ddata, cf, sfile, sformat)
          else 
            call vecio_readArray_Dble (p_Ddata, cf, sfile, sformat, p_Ipermutation)
          end if
        else
          if (.not. associated(p_Ipermutation)) then
            call vecio_readArray_Dble (p_Ddata, cf, sfile)
          else 
            call vecio_readArray_Dble (p_Ddata, cf, sfile, Ipermutation=p_Ipermutation)
          end if
        end if
      end do
    case DEFAULT
      print *,'vecio_readBlockVectorHR: Unsupported vector precision.'
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
    type(t_vectorScalar), intent(IN) :: rvector
    
    ! Name of the vector
    character(len=*), intent(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it's unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(IN) :: bunsort

    ! Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(IN) :: sformat
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
        print *, 'vecio_writeVectorHR: Could not open file '// &
                 trim(sfile)
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
      if (bunsort .and. (lsyssc_isVectorSorted (rvector))) then
        call storage_getbase_int (rvector%h_IsortPermutation,p_Ipermutation)
        ! We must use the inverse permutation
        p_Ipermutation => p_Ipermutation(rvector%NEQ+1:)
      end if

      call lsyssc_getbase_double (rvector,p_Ddata)
      
      if (.not. associated(p_Ipermutation)) then
        call vecio_writeMapleArray_Dble (p_Ddata, cf, sformat)
      else 
        call vecio_writeMapleArray_Dble (p_Ddata, cf, sformat, p_Ipermutation)
      end if
      
      ! Footer
      write (cf,'(A)') '):';
      
    case DEFAULT
      print *,'vecio_writeVectorHR: Unsupported vector precision.'
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
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! Name of the vector
    character(len=*), intent(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Write unsorted vector.
    ! =TRUE:  If the vector is sorted, it's unsorted on the fly.
    ! =FALSE: Write vector as it is.
    logical, intent(IN) :: bunsort

    ! Format string to use for the output; e.g. '(E20.10)'.
    ! If not specified, data is written to the file unformatted 
    ! (i.e. in a computer dependent, not human readable form).
    character(len=*), intent(IN) :: sformat
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
        print *, 'vecio_writeBlockVectorMaple: Could not open file '// &
                 trim(sfile)
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
        if (bunsort .and. (lsyssc_isVectorSorted (rvector%RvectorBlock(iblock)))) then
          call storage_getbase_int (&
            rvector%RvectorBlock(iblock)%h_IsortPermutation,p_Ipermutation)
          ! We must use the inverse permutation
          p_Ipermutation => p_Ipermutation(rvector%RvectorBlock(iblock)%NEQ+1:)
        end if

        call lsyssc_getbase_double (rvector%RvectorBlock(iblock),p_Ddata)
        
        if (.not. associated(p_Ipermutation)) then
          call vecio_writeMapleArray_Dble (p_Ddata, cf, sformat)
        else 
          call vecio_writeMapleArray_Dble (p_Ddata, cf, sformat, p_Ipermutation)
        end if
        
        ! If this is not the last block, attach more data
        if (iblock .ne. rvector%nblocks) then
          write (cf,'(A)', ADVANCE='NO') ','
        end if
        
      end do
      ! Footer
      write (cf,'(A)') '):';
      
    case DEFAULT
      print *,'vecio_writeBlockVectorMaple: Unsupported vector precision.'
      call sys_halt()
    end select
    
    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)

  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine vecio_writeMapleArray_Dble (Ddata, ifile, sformat, Ipermutation)

!<description>  
  ! INTERNAL SUBROUTINE.
  ! Writes the data of an array to the Maple output file iodentified by the
  ! output channel ifile.
!</description>
  
!<input>
  ! vector: array [:] of double
  real(DP), dimension(:), intent(IN) :: Ddata
  
  ! output channel to use for output
  integer, intent(IN) :: ifile
  
  ! Format string to use for the output; e.g. '(E20.10)'.
  ! If not specified, data is written to the file unformatted 
  ! (i.e. in a computer dependent, not human readable form).
  character(len=*), intent(IN) :: sformat
  
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

end module
