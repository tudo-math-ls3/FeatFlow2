!##############################################################################
!# ****************************************************************************
!# <name> boundarycondaux </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!#
!# This module provides auxiliary data structures and generic routines
!# for the treatment of boundary conditions.
!#
!#
!# The following routines are available:
!#
!# 1.) bdrc_readBoundaryCondition
!#     -> Reads boundary conditions from parameter file
!#
!# 2.) bdrc_release
!#     -> Releases a set of boundary conditions
!#
!# 3.) bdrc_calcMatrixPeriodic = bdrc_calcMatrixScalarPeriodic
!#     -> Calculates the matrix for periodic boundary conditions
!#
!# 4.) bdrc_infoBoundaryCondition
!#      -> Outputs information about the boundary condition (mostly used for debugging)
!#
!# 5.) bdrc_getNearestNeighbor2d
!#     -> Calculates number of the boundary vertex which is the
!#        nearest neighbor to some parameter value
!#
!# 6.) bdrc_createRegion
!#     -> Get the characteristics of a boundary segment and create
!#        a boundary region structure from it.
!#
!# 7.) bdrc_getNumberOfRegions
!#     -> Get the number of boundary regions on the boundary or some component
!# </purpose>
!##############################################################################

module boundarycondaux

#include "flagship.h"

!$use omp_lib
  use basicgeometry
  use boundary
  use fparser
  use fsystem
  use genoutput
  use io
  use linearsystemblock
  use linearsystemscalar
  use matrixmodification
  use storage
  use triangulation

  implicit none

  private
  public :: t_boundaryCondition
  public :: bdrc_calcMatrixPeriodic
  public :: bdrc_readBoundaryCondition
  public :: bdrc_release
  public :: bdrc_infoBoundaryCondition
  public :: bdrc_createRegion
  public :: bdrc_getNearestNeighbor2d
  public :: bdrc_getNumberOfRegions

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface bdrc_calcMatrixPeriodic
    module procedure bdrc_calcMatrixScalarPeriodic
  end interface

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>

!<constantblock>

  ! Maximum length of a section name.
  integer, parameter, public :: BDRC_MLSECTION = 64
!</constantblock>

!<constantblock description="Flags for boundary condition specification bitfield">

  ! Impose boundary condition in strong sense by filtering the system
  ! matrix and the right-hand side/residul vector explicitly
  integer(I32), parameter, public :: BDRC_STRONG   = 2**8

  ! Impose boundary condition in weak sense by overwriting the
  ! integral terms resulting from integration by parts and introducing
  ! additional integral terms into the bilinear and linear forms
  integer(I32), parameter, public :: BDRC_WEAK     = 2**9

  ! Adopt lumped assembly of the bilinear form
  integer(I32), parameter, public :: BDRC_LUMPED   = 2**10

  ! Bitmask for boundary condition (bits 0..7): Each application can
  ! defined up to 255 different tpyes of boundary conditions which can
  ! be addressed by the numbers 1..255
  integer(I32), parameter, public :: BDRC_TYPEMASK = 255
!</constantblock>

!<constantblock description="Types of special boundary conditions">

  ! Periodic boundary condition (symmetric)
  ! This condition couples two boundary segments periodically
  integer, parameter, public :: BDRC_PERIODIC      = 254

  ! Periodic boundary condition (anti-symmetric)
  ! This condition couples two boundary segments periodically
  integer, parameter, public :: BDRC_ANTIPERIODIC  = 255
!</constantblock>
!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>

!<typeblock>

  ! This data structure stores the continuous boundary conditions

  type t_boundaryCondition

    ! Number of boundary components
    integer :: iboundarycount = -1

    ! Number of spatial dimensions
    integer :: ndimension = 0

    ! Maximum number of boundary expressions
    integer :: nmaxExpressions = -1

    ! Specifier for periodic boundary conditions
    logical :: bPeriodic = .false.

    ! Specifier for strong boundary conditions
    logical :: bStrongBdrCond = .false.

    ! Specifier for weak boundary conditions
    logical :: bWeakBdrCond = .false.

    ! Handle to
    !     p_IbdrCondCpIdx = array [1..NBCT+1]
    ! which stores boundary component index vector.
    integer :: h_IbdrCondCpIdx = ST_NOHANDLE

    ! Handle to
    !     p_DmaxPar = array [1..NNCOMP]
    ! which stores the maximum parameter value of each boundary segment
    integer :: h_DmaxPar = ST_NOHANDLE

    ! Handle to
    !     p_IsegCount = array [1..NBCT]
    ! which stores the number of segments per boundary component
    integer :: h_IsegCount = ST_NOHANDLE
    
    ! Handle to
    !     p_IbdrCondType = array [1..NNCOMP]
    ! which stores the type of boundary condition of each boundary segment
    integer :: h_IbdrCondType = ST_NOHANDLE

    ! Handle to
    !     p_BisSegClosed = array [1..NNCOMP]
    ! which is .TRUE. if the right endpoint of the boundary segment also
    ! belongs to the boundary segment and .FALSE. otherwise.
    integer :: h_BisSegClosed = ST_NOHANDLE

    ! Handle to
    !    p_IbdrCompPeriodic = array [1..NNCOMP]
    ! which stores the component number of the periodic boundary neighbor
    integer :: h_IbdrCompPeriodic = ST_NOHANDLE

    ! Handle to
    !    p_IbdrCondPeriodic = array [1..NNCOMP]
    ! which stores the segment number of the periodic boundary neighbor
    integer :: h_IbdrCondPeriodic = ST_NOHANDLE

    ! Function parser for evaluation of boundary values
    type(t_fparser) :: rfparser

  end type t_boundaryCondition

!</typeblock>
!</types>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

contains

!<subroutine>

  subroutine bdrc_readBoundaryCondition(rboundaryCondition, sfilename,&
      ssectionname, ndimension, fcb_parseBoundaryCondition)

!<description>
    ! This subroutine reads boundary conditions from the parameter
    ! file and generates all internal data structures.
!</description>

!<input>
    ! The name of the parameter file to read
    character(LEN=*), intent(in) :: sfilename

    ! The name of the section in the parameter file
    character(LEN=*), intent(in) :: ssectionname

    ! The number of spatial dimensions
    integer, intent(in) :: ndimension

    ! Callback function to compute the number of expressions
    ! required for a particular type of boundary condition
    include 'intf_bdrcallback.inc'
!</input>

!<output>
    ! The boundary conditions
    type(t_boundarycondition), intent(out) :: rboundarycondition
!</output>
!</subroutine>

    ! Local variables
    real(DP), dimension(:), pointer :: p_DmaxPar
    integer, dimension(:), pointer :: p_IsegCount
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic
    integer, dimension(:), pointer :: p_IbdrCondPeriodic
    logical, dimension(:), pointer :: p_BisSegClosed
    character(FPAR_VARLEN), dimension(:), allocatable :: Svariables
    character(FPAR_STRLEN), dimension(:), allocatable :: Sexpressions
    character(FPAR_STRLEN) :: skeyword,sdata,svalue,svariable
    character(BDRC_MLSECTION) :: sectionname
    integer :: ibct,ibct1,icomp,idatalen,iexpr,ios,ipos,iunit,ivar,jpos,kpos,nexpr,nncomp

    ! Set spatial dimension
    rboundaryCondition%ndimension = ndimension

    ! Open the file
    call io_openFileForReading(sfilename, iunit, .true.)

    ! Oops...
    if (iunit .eq. -1) then
      call output_line('Unable to open input file!',&
          OU_CLASS_WARNING,OU_MODE_STD,'bdrc_readBoundaryCondition')
      call sys_halt()
    end if
    
    sectionname = adjustl(ssectionname)
    call sys_toupper(sectionname)

    ! Read through the complete file and look for the specified
    ! boundary section section
    ios = 0
    do while(ios .eq. 0)
      
      ! Read next line in file
      call io_readlinefromfile(iunit, sdata, idatalen, ios)

      ! Check for opening brackets
      ipos = scan(sdata(1:idatalen), "[")
      if (ipos .eq. 0) cycle

      ! Check for closing brackets
      jpos = scan(sdata(ipos:idatalen), "]")
      if (jpos .le. ipos) cycle

      ! Check for "[SECTIONNAME]"
      call sys_toupper(sdata(ipos+1:jpos-1), skeyword)
      if (trim(adjustl(skeyword)) .eq. trim(adjustl(sectionName))) then

        !-----------------------------------------------------------------------
        ! Read in boundary conditions from section [SECTIONNAME]
        !-----------------------------------------------------------------------

        ! (1) Read number of boundary components 'NBCT'
        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
              OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if

        call sys_toupper(sdata(1:idatalen), skeyword)        
        if (trim(adjustl(skeyword)) .ne. 'NBCT') then
          call output_line('NBCT missing!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if

        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
              OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if
        
        ! Store value of NBCT in boundary condition structure
        read(sdata(1:idatalen),*) rboundaryCondition%iboundarycount
        
        ! (2) Read maximum number of boundary expressions 'NEXPR'
        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
              OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if

        call sys_toupper(sdata(1:idatalen), skeyword)        
        if (trim(adjustl(skeyword)) .ne. 'NEXPR') then
          call output_line('NEXPR missing!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if

        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
              OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if
        
        ! Store value of NEXPR in boundary condition structure
        read(sdata(1:idatalen),*) rboundaryCondition%nmaxExpressions
        
        ! Allocate an array containing the pointers to boundary components
        call storage_new('bdrc_readBoundaryCondition', 'h_IbdrCondCpIdx',&
            rboundaryCondition%iboundarycount+1, ST_INT,&
            rboundaryCondition%h_IbdrCondCpIdx, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
        
        ! Allocate an array containing the numbers of boundary
        ! segments at each boundary components
        call storage_new('bdrc_readBoundaryCondition', 'h_IsegCount',&
            rboundaryCondition%iboundarycount, ST_INT,&
            rboundaryCondition%h_IsegCount, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(rboundaryCondition%h_IsegCount, p_IsegCount)
        
        ! Initialise the number of components
        nncomp = 0
        
        ! Loop over all boundary components
        do ibct = 1, rboundaryCondition%iboundarycount
          
          ! Set index for first boundary segment of component IBCT
          p_IbdrCondCpIdx(ibct) = nncomp+1

          ! (3) Read 'IBCT'
          call io_readlinefromfile(iunit, sdata, idatalen, ios)
          if (ios .ne. 0) then
            call output_line('Syntax error in input file!',&
                OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if
          
          call sys_toupper(sdata(1:idatalen), skeyword)        
          if (trim(adjustl(skeyword)) .ne. 'IBCT') then
            call output_line('IBCT missing!',&
                OU_CLASS_ERROR,OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if
          
          call io_readlinefromfile(iunit, sdata, idatalen, ios)
          if (ios .ne. 0) then
            call output_line('Syntax error in input file!',&
                OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if
          
          ! Compare read IBCT with current IBCT
          read(sdata(1:idatalen),*) ibct1
          if (ibct .ne. ibct1) then
            call output_line('Conflict with IBCT while reading boundary conditions!',&
                OU_CLASS_ERROR,OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if

          ! (4) Read 'NCOMP'
          call io_readlinefromfile(iunit, sdata, idatalen, ios)
          if (ios .ne. 0) then
            call output_line('Syntax error in input file!',&
                OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if
          
          call sys_toupper(sdata(1:idatalen), skeyword)        
          if (trim(adjustl(skeyword)) .ne. 'NCOMP') then
            call output_line('NCOMP missing!',&
                OU_CLASS_ERROR,OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if

          call io_readlinefromfile(iunit, sdata, idatalen, ios)
          if (ios .ne. 0) then
            call output_line('Syntax error in input file!',&
                OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if

          ! Read NCOMP and increment component counter
          read(sdata(1:idatalen),*) p_IsegCount(ibct)
          nncomp = nncomp+p_IsegCount(ibct) 
        end do ! ibct

        ! Set index of last boundary segment
        p_IbdrCondCpIdx(rboundaryCondition%iboundarycount+1) = nncomp+1
        
        ! Allocate data arrays
        call storage_new('bdrc_readBoundaryCondition', 'h_DmaxPar',&
            nncomp, ST_DOUBLE, rboundaryCondition%h_DmaxPar,&
            ST_NEWBLOCK_NOINIT)
        call storage_new('bdrc_readBoundaryCondition', 'h_IbdrCondType',&
            nncomp, ST_INT, rboundaryCondition%h_IbdrCondType,&
            ST_NEWBLOCK_NOINIT)
        call storage_new('bdrc_readBoundaryCondition', 'h_IbdrCompPeriodic',&
            nncomp, ST_INT, rboundaryCondition%h_IbdrCompPeriodic,&
            ST_NEWBLOCK_ZERO)
        call storage_new('bdrc_readBoundaryCondition', 'h_IbdrCondPeriodic',&
            nncomp, ST_INT, rboundaryCondition%h_IbdrCondPeriodic,&
            ST_NEWBLOCK_ZERO)
        call storage_new('bdrc_readBoundaryCondition', 'h_BisSegClosed',&
            nncomp, ST_LOGICAL, rboundaryCondition%h_BisSegClosed,&
            ST_NEWBLOCK_NOINIT)
        
        ! Set pointers
        call storage_getbase_double(rboundaryCondition%h_DmaxPar,&
            p_DmaxPar)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondType,&
            p_IbdrCondType)
        call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic,&
            p_IbdrCompPeriodic)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic,&
            p_IbdrCondPeriodic)
        call storage_getbase_logical(rboundaryCondition%h_BisSegClosed,&
            p_BisSegClosed)
        
        ! Initialise parser for mathematical expressions
        call fparser_create(rboundaryCondition%rfparser,&
            nncomp*rboundaryCondition%nmaxExpressions)
        
        ! (5) Read boundary parameter intervals, type of boundary and
        !     boundary expression from parameter file
        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
              OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if
        
        call sys_toupper(sdata(1:idatalen), skeyword)        
        if (trim(adjustl(skeyword)) .ne. 'PARAMETERS') then
          call output_line('PARAMETERS missing!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrc_readBoundaryCondition')
          call sys_halt()
        end if

        ! Loop over all components
        do icomp = 1, nncomp

          ! Read next line from file
          call io_readlinefromfile(iunit, sdata, idatalen, ios)
          if (ios .ne. 0) then
            call output_line('Syntax error in input file!',&
                OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if
          
          ! Split current line into individual tokens
          read(sdata(1:idatalen),*) p_DmaxPar(icomp), p_BisSegClosed(icomp), skeyword

          ! Parse type of boundary conditions and get the number of
          ! mathematical expressions required
          call fcb_parseBoundaryCondition(skeyword, ndimension,&
              p_IbdrCondType(icomp), nexpr)

          ! Check if the number of boundary conditions exceeds the
          ! maximum number of expressions
          if (nexpr .gt. rboundaryCondition%nmaxExpressions) then
            call output_line('Number of boundary conditions exceeds NEXPR!',&
                OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
            call sys_halt()
          end if

          ! Set indicator for strong boundary conditions
          if (iand(int(p_IbdrCondType(icomp),I32), BDRC_STRONG) .eq. BDRC_STRONG)&
              rboundaryCondition%bStrongBdrCond = .true.
          
          ! Set indicator for weak boundary conditions
          if (iand(int(p_IbdrCondType(icomp),I32), BDRC_WEAK) .eq. BDRC_WEAK)&
              rboundaryCondition%bWeakBdrCond = .true.
          
          ! How many mathematical expressions are required for
          ! this type of boundary conditions?
          
          ! Determine starting position of boundary data
          do ipos = 1, idatalen-len_trim(skeyword)
            if (sdata(ipos:ipos+len_trim(skeyword)-1)&
                .eq. trim(skeyword)) exit
          end do
          
          ! Copy remainder of sdata to svalue
          if (sdata(ipos+len_trim(skeyword):ipos+len_trim(skeyword)) .eq. "'"  .or.&
              sdata(ipos+len_trim(skeyword):ipos+len_trim(skeyword)) .eq. '"') then
            svalue = trim(adjustl(sdata(ipos+len_trim(skeyword)+1:)))
          else
            svalue = trim(adjustl(sdata(ipos+len_trim(skeyword):)))
          end if
          
          ! Concatenate multi-line expressions
          do while(ios .eq. 0)
            ! Get length of expression
            ipos = len_trim(svalue)
            
            ! Check if expression is continued in the following line
            if (svalue(max(1, ipos-2):ipos) .eq. '...') then
              ipos = ipos-2
            else
              exit
            end if
            
            ! Read next line in file
            call io_readlinefromfile(iunit, sdata, idatalen, ios)
            if (ios .ne. 0) then
              call output_line('Syntax error in input file!',&
                  OU_CLASS_ERROR, OU_MODE_STD,'bdrc_readBoundaryCondition')
              call sys_halt()
            end if
            
            ! Append line
            svalue(ipos:) = trim(adjustl(sdata(1:idatalen)))
          end do
          
          ! Extract the symbolic variables (if any)
          jpos = scan(svalue(1:ipos), ";", .true.)
          
          if (jpos .eq. 0) then
            
            ! Allocate temporal memory
            allocate(Sexpressions(min(rboundaryCondition%nmaxExpressions,nexpr)))
            
            ! Split mathematical expressions into tokens
            read(svalue(1:ipos),*) Sexpressions
            
            ! Loop over all expressions and apply them to function
            ! parser, whereby no symbolic variables are given
            do iexpr = 1, min(nexpr,rboundaryCondition%nmaxExpressions)
              call fparser_parseFunction(rboundaryCondition%rfparser,&
                  rboundaryCondition%nmaxExpressions*(icomp-1)+iexpr,&
                  Sexpressions(iexpr), (/''/))
            end do
            
            ! Deallocate temporal memory
            deallocate(Sexpressions)
            
          else
            
            ! Counte number of symbolic variables
            svariable = svalue(jpos+1:ipos)
            
            kpos = scan(svariable, ","); ivar = 0
            do while (kpos .ne. 0)
              ivar = ivar+1
              svariable = trim(adjustl(svariable(kpos+1:len_trim(svariable))))
              kpos = scan(svariable, ",")
            end do
            
            ! Allocate temporal memory
            allocate(Svariables(ivar+1),&
                     Sexpressions(min(nexpr,rboundaryCondition%nmaxExpressions)))
            
            ! Split mathematical expressions into tokens
            read(svalue(1:ipos),*) Sexpressions
            
            ! Initialise symbolic variables
            svariable = svalue(jpos+1:ipos)
            
            kpos = scan(svariable, ","); ivar = 0
            do while (kpos .ne. 0)
              ivar = ivar+1
              Svariables(ivar) = trim(adjustl(svariable(1:kpos-1)))
              svariable = trim(adjustl(svariable(kpos+1:len_trim(svariable))))
              kpos = scan(svariable, ",")
            end do
            Svariables(ivar+1) = trim(adjustl(svariable(1:len_trim(svariable))))
            
            ! Loop over all expressions and apply them to function parser
            do iexpr = 1, min(nexpr,rboundaryCondition%nmaxExpressions)
              call fparser_parseFunction(rboundaryCondition%rfparser,&
                  rboundaryCondition%nmaxExpressions*(icomp-1)+iexpr,&
                  Sexpressions(iexpr), Svariables)
            end do
            
            ! Deallocate temporal memory
            deallocate(Svariables,Sexpressions)
          end if

          ! Initialise empty expressions by zero
          do iexpr = max(1,nexpr+1), rboundaryCondition%nmaxExpressions
            call fparser_parseFunction(rboundaryCondition%rfparser,&
                rboundaryCondition%nmaxExpressions*(icomp-1)+iexpr,&
                '0', (/''/))
          end do

          ! Do we have periodic boundary conditions
          if ((iand(int(p_IbdrCondType(icomp),I32),&
                    BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
              (iand(int(p_IbdrCondType(icomp),I32),&
                    BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
            
            ! Allocate temporal memory
            allocate(Sexpressions(nexpr+2))
            read(svalue(1:ipos),*) Sexpressions
            
            ! The two last tokens specify the component and segment of
            ! the periodic boundary part, respectively.
            read(Sexpressions(nexpr+1),*) p_IbdrCompPeriodic(icomp)
            read(Sexpressions(nexpr+2),*) p_IbdrCondPeriodic(icomp)

            ! Mark presence of periodic boundary conditions
            rboundaryCondition%bPeriodic = .true.
            
            ! Deallocate temporal memory
            deallocate(Sexpressions)
          end if

        end do ! icomp

        ! That is all
        close(iunit)
        return
      end if
    end do
    
    ! Close the file, finish
    close(iunit)

    call output_line('Unable to read boundary conditions from file!',&
        OU_CLASS_ERROR,OU_MODE_STD,'bdrc_readBoundaryCondition')
    call sys_halt()
    
  end subroutine bdrc_readBoundaryCondition

  ! *****************************************************************************

!<subroutine>

  subroutine bdrc_release(rboundaryCondition)

!<description>
    ! This subroutine releases a boundary condition
!</description>

!<inputoutput>
    ! The boundary conditions
    type(t_boundaryCondition), intent(inout) :: rboundaryCondition
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rboundaryCondition%h_IbdrCondCpIdx .ne. ST_NOHANDLE)&
        call storage_free(rboundaryCondition%h_IbdrCondCpIdx)
    if (rboundaryCondition%h_DmaxPar .ne. ST_NOHANDLE)&
        call storage_free(rboundaryCondition%h_DmaxPar)
    if (rboundaryCondition%h_IsegCount .ne. ST_NOHANDLE)&
        call storage_free(rboundaryCondition%h_IsegCount)
    if (rboundaryCondition%h_IbdrCondType .ne. ST_NOHANDLE)&
        call storage_free(rboundaryCondition%h_IbdrCondType)
    if (rboundaryCondition%h_IbdrCompPeriodic .ne. ST_NOHANDLE)&
        call storage_free(rboundaryCondition%h_IbdrCompPeriodic)
    if (rboundaryCondition%h_IbdrCondPeriodic .ne. ST_NOHANDLE)&
        call storage_free(rboundaryCondition%h_IbdrCondPeriodic)
    if (rboundaryCondition%h_BisSegClosed .ne. ST_NOHANDLE)&
        call storage_free(rboundaryCondition%h_BisSegClosed)

    ! Release function parser
    call fparser_release(rboundaryCondition%rfparser)

  end subroutine bdrc_release

  ! *****************************************************************************

!<subroutine>

  subroutine bdrc_calcMatrixScalarPeriodic(rboundaryCondition,&
      rmatrix, rtriangulation)

!<description>
    ! This subroutine calculates the modified matrix structure
    ! required to impose periodic boundary conditions in 1D and 2D.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation
!</input>

!<inputoutput>
    ! Scalar matrix to be adjusted
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:), pointer :: p_DmaxPar
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:,:), pointer :: p_Irows
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic
    integer, dimension(:), pointer :: p_IbdrCondPeriodic
    logical, dimension(:), pointer :: p_BisSegClosed
    integer, dimension(2) :: Isize
    integer :: h_Irows,nrows


    ! Check if periodic boundary conditions are present
    if (rboundaryCondition%bPeriodic) then

      ! Create memory for pairs of boundary vertices
      Isize = (/2,rtriangulation%NVBD/)
      call storage_new('bdrc_calcMatrixScalarPeriodic2D', 'p_Irows',&
          Isize, ST_INT, h_Irows, ST_NEWBLOCK_ZERO)
      call storage_getbase_int2D(h_Irows, p_Irows)

      ! Get underlying triangulation structure
      if (present(rtriangulation)) then
        p_rtriangulation => rtriangulation
      else
        if (.not.associated(rmatrix%p_rspatialDiscrTrial)) then
          call output_line('No discretisation associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrc_calcMatrixScalarPeriodic')
          call sys_halt()
        end if

        if (.not.associated(rmatrix%p_rspatialDiscrTrial%p_rtriangulation)) then
          call output_line('No triangulation associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrc_calcMatrixScalarPeriodic')
          call sys_halt()
        end if
        p_rtriangulation => rmatrix%p_rspatialDiscrTrial%p_rtriangulation
      end if

      ! Check spatial dimensions
      if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
        call output_line('Spatial dimension mismatch!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrc_calcMatrixScalarPeriodic')
        call sys_halt()
      end if

      ! How many spatial dimensions do we have?
      select case (rboundaryCondition%ndimension)
      case (NDIM1D)
        ! Set pointers for triangulation
        call storage_getbase_int(p_rtriangulation%h_IboundaryCpIdx,&
            p_IboundaryCpIdx)
        call storage_getbase_int(p_rtriangulation&
            %h_IverticesAtBoundary, p_IverticesAtBoundary)

        ! Set pointers for boundary
        call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
        call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)

        ! Calculate pairs of vertices in 1D
        call calcPeriodic_1D(p_IbdrCompPeriodic, p_IbdrCondPeriodic,&
            p_IbdrCondType, p_IbdrCondCpIdx, rboundaryCondition&
            %iboundarycount, p_IverticesAtBoundary, p_IboundaryCpIdx, p_Irows, nrows)

      case (NDIM2D)
        ! Set pointers for triangulation
        call storage_getbase_double(&
            p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
        call storage_getbase_int(&
            p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
        call storage_getbase_int(&
            p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

        ! Set pointers for boundary
        call storage_getbase_double (rboundaryCondition%h_DmaxPar, p_DmaxPar)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
        call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
        call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)
        call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

        ! Calculate pairs of vertices in 2D
        call calcPeriodic_2D(p_IbdrCompPeriodic, p_IbdrCondPeriodic,&
            p_IbdrCondType, p_IbdrCondCpIdx, p_DmaxPar,&
            p_BisSegClosed, rboundaryCondition%iboundarycount,&
            p_IverticesAtBoundary, p_IboundaryCpIdx,&
            p_DvertexParameterValue, p_Irows, nrows)

      case default
        call output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrc_calcMatrixScalarPeriodic')
        call sys_halt()
      end select

      ! Modify the matrix accordingly and require symmetric sparsity
      call mmod_mergeLines(rmatrix, p_Irows(:,1:nrows), .true.)

      ! Clear temporal memory
      call storage_free(h_Irows)
    end if

  contains

    ! Here are the real working routines

    !***************************************************************
    ! Calculates the pairs of periodic boundary vertices in 1D
    subroutine calcPeriodic_1D(IbdrCompPeriodic, IbdrCondPeriodic,&
        IbdrCondType, IbdrCondCpIdx, nbct, IverticesAtBoundary,&
        IboundaryCpIdx, Irows, nrows)

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Auxiliary array for pairs of boundary vertices
      integer, dimension(:,:), intent(out) :: Irows

      ! Number of pairs of boundary vertices
      integer, intent(out) :: nrows

      ! local variables
      integer :: ivbd,ivbdPeriodic,ibct,ibctPeriodic,isegment


      ! Initialise row counter
      nrows = 0

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Are we strongly imposed periodic boundary conditions?
        if (iand(int(IbdrCondType(isegment),I32), BDRC_STRONG)   .eq. BDRC_STRONG .and.&
           (iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK) .eq. BDRC_PERIODIC .or.&
            iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then

          ! Compute vertex parameter value at periodic boundary
          ibctPeriodic = IbdrCompPeriodic(isegment)
          ivbdPeriodic = IboundaryCpIdx(ibctPeriodic)
          
          ! Append pair of periodic vertices
          nrows = nrows+1
          Irows(1, nrows) = IverticesAtBoundary(ivbd)
          Irows(2, nrows) = IverticesAtBoundary(ivbdPeriodic)
        end if

      end do
    end subroutine calcPeriodic_1D


    !***************************************************************
    ! Calculates the pairs of periodic boundary vertices in 2D
    subroutine calcPeriodic_2D(IbdrCompPeriodic, IbdrCondPeriodic,&
        IbdrCondType, IbdrCondCpIdx, DmaxParam, BisSegClosed, nbct,&
        IverticesAtBoundary, IboundaryCpIdx, DvertexParameterValue,&
        Irows, nrows)

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Auxiliary array for pairs of boundary vertices
      integer, dimension(:,:), intent(out) :: Irows

      ! Number of pairs of boundary vertices
      integer, intent(out) :: nrows

      ! local variables
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic
      integer :: ibct,ibctPeriodic,isegment,isegmentPeriodic


      ! Initialise row counter
      nrows = 0

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue =&
            nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1) isegment =&
                IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Are we strongly imposed periodic boundary conditions?
          if (iand(int(IbdrCondType(isegment),I32), BDRC_STRONG)   .eq. BDRC_STRONG .and.&
              iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then

            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            if (isegmentPeriodic .eq. IbdrCondCpIdx(ibctPeriodic)) then
              dVertexParameterPeriodic = DmaxParam(isegment)-&
                  DvertexParameterValue(ivbd)
            else
              dVertexParameterPeriodic = DmaxParam(isegmentPeriodic-1)+&
                  DmaxParam(isegment)-DvertexParameterValue(ivbd)
            end if

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast,&
                ivbdPeriodic)

            ! Append pair of periodic vertices
            nrows = nrows+1
            Irows(1, nrows) = IverticesAtBoundary(ivbd)
            Irows(2, nrows) = IverticesAtBoundary(ivbdPeriodic)

          elseif (iand(int(IbdrCondType(isegment),I32), BDRC_STRONG)   .eq. BDRC_STRONG .and.&
                  iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC) then

            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            dVertexParameterPeriodic = DmaxParam(isegmentPeriodic)-&
                (DmaxParam(isegment)-DvertexParameterValue(ivbd))

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast,&
                ivbdPeriodic)

            ! Append pair of periodic vertices
            nrows = nrows+1
            Irows(1, nrows) = IverticesAtBoundary(ivbd)
            Irows(2, nrows) = IverticesAtBoundary(ivbdPeriodic)

          end if

        end do
      end do
    end subroutine calcPeriodic_2D
  end subroutine bdrc_calcMatrixScalarPeriodic

  ! *****************************************************************************

!<subroutine>

  subroutine bdrc_infoBoundaryCondition(rboundaryCondition)

!<description>
    ! This subroutine output information about the boundary condition
!</description>

!<input>
    ! The boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition
!</input>
!</subroutine

    ! local variables
    real(DP), dimension(:), pointer :: p_DmaxPar
    integer, dimension(:), pointer :: p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCompPeriodic
    integer, dimension(:), pointer :: p_IbdrCondPeriodic
    logical, dimension(:), pointer :: p_BisSegClosed
    integer :: ibct,icomp


    call output_line ('BoundaryCondition:')
    call output_line ('------------------')
    call output_line ('iboundarycount:     '//&
        trim(sys_siL(rboundaryCondition%iboundarycount,15)))
    call output_line ('ndimension:         '//&
        trim(sys_siL(rboundaryCondition%ndimension,15)))
    call output_line ('nmaxExpressions:    '//&
        trim(sys_siL(rboundaryCondition%nmaxExpressions,15)))
    call output_line ('bPeriodic:          '//&
        merge('true ','false',rboundaryCondition%bPeriodic))
    call output_line ('bStrongBdrCond:     '//&
        merge('true ','false',rboundaryCondition%bStrongBdrCond))
    call output_line ('bWeakBdrCond:       '//&
        merge('true ','false',rboundaryCondition%bWeakBdrCond))
    call output_line ('h_IbdrCondCpIdx:    '//&
        trim(sys_siL(rboundaryCondition%h_IbdrCondCpIdx,15)))
    call output_line ('h_DmaxPar:          '//&
        trim(sys_siL(rboundaryCondition%h_DmaxPar,15)))
    call output_line ('h_IbdrCondType:     '//&
        trim(sys_siL(rboundaryCondition%h_IbdrCondType,15)))
    call output_line ('h_BisSegClosed:     '//&
        trim(sys_siL(rboundaryCondition%h_BisSegClosed,15)))
    call output_line ('h_IbdrCompPeriodic: '//&
        trim(sys_siL(rboundaryCondition%h_IbdrCompPeriodic,15)))
    call output_line ('h_IbdrCondPeriodic: '//&
        trim(sys_siL(rboundaryCondition%h_IbdrCondPeriodic,15)))
    call output_lbrk()

    ! Set pointers
    call storage_getbase_double(rboundaryCondition%h_DmaxPar,&
        p_DmaxPar)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondType,&
        p_IbdrCondType)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx,&
        p_IbdrCondCpIdx)
    call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic,&
        p_IbdrCompPeriodic)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic,&
        p_IbdrCondPeriodic)
    call storage_getbase_logical(rboundaryCondition%h_BisSegClosed,&
        p_BisSegClosed)

    ! Loop over boundary components
    do ibct = 1, rboundaryCondition%iboundarycount

      call output_lbrk()
      call output_separator(OU_SEP_PLUS)
      call output_line ('Boundary component: '//trim(sys_siL(ibct,15)))
      call output_separator(OU_SEP_PLUS)
      call output_lbrk()
      
      ! Loop over all boundary segments in current boundary component
      do icomp = p_IbdrCondCpIdx(ibct), p_IbdrCondCpIdx(ibct+1)-1

        if (icomp .eq. p_IbdrCondCpIdx(ibct)) then
          
          ! First segment needs some special treatment
          call output_line ('Segment: '//trim(sys_siL(icomp,15)))
          call output_separator(OU_SEP_MINUS)
          call output_line('Type of boundary condition:                '//&
              trim(sys_siL(p_IbdrCondType(icomp),15)))
          if (rboundaryCondition%ndimension .eq. 2) then
            if (p_BisSegClosed(icomp)) then
              call output_line('Type of boundary segment in 2D:            '//&
                  merge('( ]','[ ]', p_BisSegClosed(p_IbdrCondCpIdx(ibct+1)-1)))
            else
              call output_line('Type of boundary segment in 2D:            '//&
                  merge('( )','[ )', p_BisSegClosed(p_IbdrCondCpIdx(ibct+1)-1)))
            end if
          end if
          call output_line('Component number of the periodic neighbor: '//&
              trim(sys_siL(p_IbdrCompPeriodic(icomp),15)))
          call output_line('Segment number of the periodic neighbor:   '//&
              trim(sys_siL(p_IbdrCondPeriodic(icomp),15)))
          call output_separator(OU_SEP_MINUS)
        else
          
          ! All other segments have a predecessor
          call output_line ('Segment: '//trim(sys_siL(icomp,15)))
          call output_separator(OU_SEP_MINUS)
          call output_line('Type of boundary condition:                '//&
              trim(sys_siL(p_IbdrCondType(icomp),15)))
          if (rboundaryCondition%ndimension .eq. 2) then
            if (p_BisSegClosed(icomp)) then
              call output_line('Type of boundary segment in 2D:            '//&
                  merge('( ]','[ ]', p_BisSegClosed(icomp-1)))
            else
              call output_line('Type of boundary segment in 2D:            '//&
                  merge('( )','[ )', p_BisSegClosed(icomp-1)))
            end if
          end if
          call output_line('Component number of the periodic neighbor: '//&
              trim(sys_siL(p_IbdrCompPeriodic(icomp),15)))
          call output_line('Segment number of the periodic neighbor:   '//&
              trim(sys_siL(p_IbdrCondPeriodic(icomp),15)))
          call output_separator(OU_SEP_MINUS)
        end if
      end do
    end do
    
  end subroutine bdrc_infoBoundaryCondition

  ! *****************************************************************************

!<subroutine>

  pure subroutine bdrc_getNearestNeighbor2d(DvertexParameterValue,&
      dParameterValue, ivbdMin, ivbdMax, ivbd)

!<description>
    ! This subroutine determines the number of the boundary vertex which
    ! represents the nearest neighbor to the desired parameter value.
!</description>

!<input>
    ! Parameter values of all boundary vertices
    real(DP), dimension(:), intent(in) :: DvertexParameterValue

    ! Parameter value to which nearest neighbor is required
    real(DP), intent(in) :: dParameterValue

    ! Minimum number of boundary vertex
    integer, intent(in) :: ivbdMin

    ! Maximum number of boundary vertex
    integer, intent(in) :: ivbdMax
!</input>

!<inputoutput>
    ! Initial guess for the number of the boundary vertex on
    ! input. The number of the nearest neighbor on output.
    integer, intent(inout) :: ivbd
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ilow, ihigh, imid

    ilow  = ivbdMin
    ihigh = ivbdMax

    if (ivbd .ge. ilow .and. ivbd .le. ihigh) then
      imid  = ivbd
    else
      imid = (ilow+ihigh)/2
    end if

    do while( ilow .le. ihigh)
      if (DvertexParameterValue(imid) .gt. dParameterValue) then
        ihigh = imid-1
      elseif (DvertexParameterValue(imid) .lt. dParameterValue) then
        ilow = imid+1
      else
        exit
      end if
      imid = (ilow+ihigh)/2
    end do
    ivbd = imid
  end subroutine bdrc_getNearestNeighbor2d

  ! *****************************************************************************

!<subroutine>

  subroutine bdrc_createRegion(rboundaryCondition, iboundCompIdx,&
      iboundSegIdx, rregion)

!<description>
    ! This subroutine creates a boundary region from a boundary
    ! condition structure which stores all required data.
!</description>

!<input>
    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! Index of boundary component.
    integer, intent(in) :: iboundCompIdx

    ! Index of the boundary segment.
    ! =0: Create a boundary region that covers the whole boundary component.
    integer, intent(in) :: iboundSegIdx
!</input>

!<output>
    ! Boundary region that is characterised by the boundary segment
    type(t_boundaryRegion), intent(out) :: rregion
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DmaxPar
    integer, dimension(:), pointer :: p_IsegCount
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    logical, dimension(:), pointer :: p_BisSegClosed
    real(DP) :: dcurrentpar, dendpar, dmaxpar
    integer :: iproperties,isegment

    if ((iboundCompIdx .gt. rboundaryCondition%iboundarycount) .or.&
        (iboundCompIdx .lt. 0)) then
      call output_line ('iboundCompIdx out of bounds!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bdrc_createRegion')
      call sys_halt()
    endif

    ! Set pointers for boundary
    call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
    call storage_getbase_int(rboundaryCondition%h_IsegCount, p_IsegCount)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
    call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)


    if (iboundSegIdx .ne. 0) then

      if ((iboundSegIdx .gt. p_IsegCount(iboundCompIdx)) .or. (iboundSegIdx.lt.0)) then
        call output_line ('iboundSegIdx out of bounds!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bdrc_createRegion')
        call sys_halt()
      endif

      ! Get segment number in absolute values
      isegment = p_IbdrCondCpIdx(iboundCompIdx)+iboundSegIdx-1

      ! Get the start and end parameter value
      if (iboundSegIdx .eq. 1) then
        dcurrentpar = 0.0_DP
      else
        dcurrentpar = p_DmaxPar(isegment-1)
      end if
      dmaxpar = p_DmaxPar(p_IbdrCondCpIdx(iboundCompIdx+1)-1)
      dendpar = p_DmaxPar(isegment)

      ! We have an unspecified boundary segment
      rregion%isegmentType = BOUNDARY_TYPE_ANALYTIC

      ! Set interval properties
      iproperties = 0

      if (iboundSegIdx .eq. 1) then
        if (.not.p_BisSegClosed(p_IbdrCondCpIdx(iboundCompIdx+1)-1))&
            iproperties = ior(iproperties, BDR_PROP_WITHSTART)
      else
        if (.not.p_BisSegClosed(isegment-1))&
            iproperties = ior(iproperties, BDR_PROP_WITHSTART)
      end if

      if (p_BisSegClosed(isegment))&
          iproperties = ior(iproperties, BDR_PROP_WITHEND)

    else

      ! Create a boundary region that covers the whole boundary component.
      dcurrentpar = 0.0_DP
      dmaxpar     = p_DmaxPar(p_IbdrCondCpIdx(iboundCompIdx+1)-1)
      dendpar     = dmaxpar

      ! We have an unspecified boundary segment
      rregion%isegmentType = BOUNDARY_TYPE_ANALYTIC

      ! Set interval properties. The whole boundary either includes
      !  the start or the end-point.
      if (p_BisSegClosed(p_IbdrCondCpIdx(iboundCompIdx+1)-1)) then
        iproperties = BDR_PROP_WITHEND
      else
        iproperties = BDR_PROP_WITHSTART
      end if

    end if

    ! Create the boundary region structure
    rregion%cparType = BDR_PAR_01
    rregion%dminParam = dcurrentpar
    rregion%dmaxParam = dendpar
    rregion%dmaxParamBC = dmaxpar
    rregion%ctype = BDR_TP_CURVE
    rregion%iproperties = iproperties
    rregion%iboundCompIdx = iboundCompIdx
    rregion%iboundSegIdx = 0

  end subroutine bdrc_createRegion

  ! *****************************************************************************

!<function>

  function bdrc_getNumberOfRegions(rboundaryCondition, ibdc) result(nregions)

!<description>
    ! This function calculates the number of boundary regions on the
    ! boundary or on the boundary component ibdc.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! OPTIONAL: boundary component
    integer, intent(in), optional :: ibdc
!</input>

!<result>
    integer :: nregions
!</result>
!</function>

    ! local variables
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer :: ibct

    ! Set pointers
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)

    if (present(ibdc)) then
      nregions = p_IbdrCondCpIdx(ibdc+1)-p_IbdrCondCpIdx(ibdc)
    else
      nregions = 0
      do ibct = 1, rboundaryCondition%iboundarycount
        nregions = nregions + p_IbdrCondCpIdx(ibct+1)-p_IbdrCondCpIdx(ibct)
      end do
    end if

  end function bdrc_getNumberOfRegions
  
end module boundarycondaux
