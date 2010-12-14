module typedsymbol

  use fsystem
  use genoutput
  use storage
  use collection

  implicit none
  
  private
  
  ! Symbol type: Invalid symbol
  integer, parameter, public :: STYPE_INVALID = 0
  
  ! Symbol type: Integer
  integer, parameter, public :: STYPE_INTEGER = 1
  
  ! Symbol type: Double precision
  integer, parameter, public :: STYPE_DOUBLE  = 2

  ! Symbol type: string
  integer, parameter, public :: STYPE_STRING  = 3

  ! Symbol type: Variable name
  integer, parameter, public :: STYPE_VAR     = 4
  
  ! Encapsules rtyped symbols.
  type t_symbolValue
  
    ! Type of the variable. A STYPE_xxxx constant.
    integer :: ctype = STYPE_INVALID
    
    ! Integer value.
    integer :: ivalue = 0
    
    ! Double precision return value
    real(DP) :: dvalue = 0.0_DP
    
    ! String value
    integer :: ilength = 0
    character(SYS_STRLEN) :: svalue
    
    ! Name of the symbol in case it represents a variable in a collection.
    ! For variables in a collection, ctype is set to the type of the
    ! variable and Xvalue to its value; svarname is set to the variable name.
    character(SYS_NAMELEN) :: svarname = ""
    
  end type
  
  public :: t_symbolValue
  
  public :: tpsym_determineSymbolType
  public :: tpsym_parseSymbol
  public :: tpsym_saveSymbol
  public :: tpsym_setinvalid
  public :: tpsym_undefine
  public :: operator(+),operator(-),operator(*),operator(/),MOD
  public :: operator(<),operator(>),operator(<=),operator(>=)
  public :: operator(==),operator(/=),assignment(=)
  
  interface assignment(=)
    module procedure tpsym_set
    module procedure tpsym_setint
    module procedure tpsym_setdouble
    module procedure tpsym_setstring
  end interface

  interface operator(+)
    module procedure tpsym_add
  end interface

  interface operator(-)
    module procedure tpsym_subtract
  end interface

  interface operator(*)
    module procedure tpsym_multiply
  end interface

  interface operator(/)
    module procedure tpsym_divide
  end interface

  interface mod
    module procedure tpsym_divide
  end interface

  interface operator(<)
    module procedure tpsym_less
  end interface

  interface operator(>)
    module procedure tpsym_greater
  end interface
  
  interface operator(<=)
    module procedure tpsym_lesseq
  end interface

  interface operator(>=)
    module procedure tpsym_greatereq
  end interface

  interface operator(==)
    module procedure tpsym_eq
  end interface

  interface operator(/=)
    module procedure tpsym_neq
  end interface
  
contains

  ! ***************************************************************************

  !<subroutine>

  subroutine cmdprs_getSymbolSection (rcollection,ssymbol,inestlevel,ssection,bfound)

  !<description>
    ! Tries to find the name of the section containing a symbol.
  !</description>

  !<inputoutput>
    ! Collection containing symbols.
    type(t_collection), intent(inout) :: rcollection
  !</inputoutput>
  
  !<input>
    ! Name of the symbol (variable,...)
    character(len=*), intent(in) :: ssymbol
    
    ! Level of nesting
    integer, intent(in) :: inestlevel
  !</input>
  
  !<output>
    ! Name of the section.
    character(len=*), intent(out) :: ssection
    
    ! Returns if the symbol was found or not.
    logical, intent(out) :: bfound
  !</output>
  
!</subroutine>
    
    ! local variables
    integer :: ilv

    bfound = .false.
    ssection = ""
    
    ! Start searching from the deepest nesting
    do ilv = inestlevel,0,-1
      if (ilv .gt. 0) then
        ! local variables
        ssection = trim(sys_siL(ilv,10))
      else
        ! Unnamed section, global variables
        ssection = ""
      end if
      
      if (collct_queryvalue(rcollection,ssymbol,ssectionName=ssection) .gt. 0) then
        ! Found.
        bfound = .true.
        return
      end if
    end do
    
    ! Not found.
  
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  elemental subroutine tpsym_determineSymbolType (ssymbol,ctype)

  !<description>
    ! Determins the type of a symbol.
  !</description>

  !<input>
    ! Symbol to analyse
    character(len=*), intent(in) :: ssymbol
  !</input>

  !<output>
    ! Type of the symbol
    integer, intent(out) :: ctype
  !</output>
  
!</subroutine>

    if (verify(ssymbol,"+-0123456789") .eq. 0) then
      ! An integer
      ctype = STYPE_INTEGER
    else if (verify(ssymbol,"+-0123456789.E") .eq. 0) then
      ! An double precision value
      ctype = STYPE_DOUBLE
    else if (verify(sys_upcase(ssymbol),"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ_") .eq. 0) then
      if (verify(ssymbol(1:1),"0123456789") .ne. 0) then
        ! A variable name
        ctype = STYPE_VAR
      end if
    else if ((index(ssymbol,"""") .eq. 1) .and. &
             (index(ssymbol,"""",.true.) .eq. LEN_TRIM(ssymbol))) then
      ! A string
      ctype = STYPE_STRING
    else
      ! Invalid qualifier
      ctype = STYPE_INVALID
    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine tpsym_parseSymbol (ssymbol,rcollection,inestlevel,rvalue)

  !<description>
    ! Parses a string and creates a symbol from it.
  !</description>

  !<input>
    ! Symbol to analyse
    character(len=*), intent(in) :: ssymbol
    
    ! Reference to a collection for variables
    type(t_collection), intent(inout), target :: rcollection
    
    ! Nesting level of the symbol
    integer, intent(in) :: inestlevel
  !</input>

  !<output>
    ! Corresponding value.
    type(t_symbolValue), intent(out) :: rvalue
  !</output>
  
!</subroutine>
    character(len=SYS_NAMELEN) :: ssection
    logical :: bfound
    integer :: cvartype

    ! Determine the type
    call tpsym_determineSymbolType (ssymbol,rvalue%ctype)
    
    ! Parse.
    select case (rvalue%ctype)
    case (STYPE_INTEGER)
      read (ssymbol,*) rvalue%ivalue
    case (STYPE_DOUBLE )
      read (ssymbol,*) rvalue%dvalue
    case (STYPE_STRING )
      read(ssymbol,"(A)") rvalue%svalue
      rvalue%ilength = len_trim(ssymbol)
    case (STYPE_VAR)
      
      ! Initialise the variable name
      rvalue%svarname = trim(ssymbol)
      
      ! Try to get the correct type etc. from the collection.
      rvalue%ctype = STYPE_INVALID
      call cmdprs_getSymbolSection (rcollection,ssymbol,inestlevel,ssection,bfound)
      if (bfound) then
        cvartype = collct_gettype(rcollection,ssymbol,ssectionName=ssection)
        ! Fetch the value if possible
        select case (cvartype)
        case (COLLCT_INTEGER)
          rvalue%ctype = STYPE_INTEGER
          rvalue%ivalue = collct_getvalue_int (rcollection, ssymbol, ssectionName=ssection)
        case (COLLCT_REAL)
          rvalue%ctype = STYPE_DOUBLE
          rvalue%ivalue = collct_getvalue_real (rcollection, ssymbol, ssectionName=ssection)
        case (COLLCT_STRING)
          rvalue%ctype = STYPE_STRING
          call collct_getvalue_string (rcollection, ssymbol, rvalue%svalue, ssectionName=ssection)
          rvalue%ilength = len_trim(rvalue%svalue)
        end select
      end if
      
    end select
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine tpsym_saveSymbol (rvalue,inestlevel,rcollection)

  !<description>
    ! Tries to save the value corresponding to a symbol into the collection.
  !</description>

  !<input>
    ! Value to save
    type(t_symbolValue), intent(in) :: rvalue

    ! Nesting level of the symbol
    integer, intent(in) :: inestlevel
  !</input>

  !<inputoutput>
    ! Reference to a collection for variables
    type(t_collection), intent(inout), target :: rcollection
  !</inputoutput>
  
!</subroutine>
    character(len=SYS_NAMELEN) :: ssection
    logical :: bfound

    ! The symbol must be from the collection!
    if (rvalue%svarname .eq. "") return
    
    call cmdprs_getSymbolSection (rcollection,rvalue%svarname,inestlevel,ssection,bfound)
    if (bfound) then
      ! Write to the collection
      select case (rvalue%ctype)
      case (STYPE_INTEGER)
        call collct_setvalue_int (rcollection, rvalue%svarname, rvalue%ivalue, .false., ssectionName=ssection)
      case (STYPE_DOUBLE)
        call collct_setvalue_real (rcollection, rvalue%svarname, rvalue%dvalue, .false., ssectionName=ssection)
      case (STYPE_STRING)
        call collct_setvalue_string (rcollection, rvalue%svarname, rvalue%svalue, .false., ssectionName=ssection)
      end select
    end if
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine tpsym_setinvalid (rvalue)

  !<description>
    ! Declares a symbol as 'invalid'.
  !</description>

  !<inputoutput>
    ! Symbol to set.
    type(t_symbolValue), intent(inout) :: rvalue
  !</inputoutput>

!</subroutine>

    rvalue%ctype = STYPE_INVALID
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine tpsym_undefine (rvalue)

  !<description>
    ! Undefines a symbol, removes all settings from the structure;
    ! reinitialisation.
  !</description>

  !<inputoutput>
    ! Symbol to set.
    type(t_symbolValue), intent(out) :: rvalue
  !</inputoutput>

!</subroutine>

  end subroutine

  ! ***************************************************************************
  ! Operators
  ! ***************************************************************************

  ! Assignment
  subroutine tpsym_set (rdest,rsource)
  type (t_symbolValue), intent(inout) :: rdest
  type (t_symbolValue), intent(in) :: rsource
  
    ! What's the target type? If the target has a type and there is a variable
    ! name associated, we must not change the type!
    if ((rdest%ctype .eq. STYPE_INVALID) .or. (rdest%svarname .eq. "")) then
      rdest%ctype = rsource%ctype
      select case (rsource%ctype)
      case (STYPE_INTEGER)
        rdest%ivalue = rsource%ivalue
      case (STYPE_DOUBLE)
        rdest%dvalue = rsource%dvalue
      case (STYPE_STRING)
        rdest%svalue = rsource%svalue
        rdest%ilength = rsource%ilength
      end select
    else
      !hTHere must probably a type conversion be done.
      ! Use the routines below to do the assignments.
      select case (rsource%ctype)
      case (STYPE_INTEGER)
        rdest = rsource%ivalue
      case (STYPE_DOUBLE)
        rdest = rsource%dvalue
      case (STYPE_STRING)
        rdest = rsource%svalue
      case default
        ! Invalid assignment
        rdest%ctype = STYPE_INVALID
      end select
    end if
  
  end subroutine

  subroutine tpsym_setint (rdest,isource)
  type (t_symbolValue), intent(inout) :: rdest
  integer, intent(in) :: isource
  
    ! Type conversion if possible
    select case (rdest%ctype)
    case (STYPE_INTEGER) 
      rdest%ivalue = isource
    case (STYPE_DOUBLE) 
      rdest%dvalue = real(isource,DP)
    case (STYPE_INVALID) 
      ! Invalid and no variable associated? Will be an int.
      if (rdest%svarname .eq. "") then
        rdest%ctype = STYPE_INTEGER
        rdest%ivalue = isource
      else
        rdest%ctype = STYPE_INVALID
      end if
    case default
      ! Otherwise invalid value.
      rdest%ctype = STYPE_INVALID
    end select
    
  end subroutine

  subroutine tpsym_setdouble (rdest,dsource)
  type (t_symbolValue), intent(inout) :: rdest
  real(DP), intent(in) :: dsource
  
    ! Type conversion if possible
    select case (rdest%ctype)
    case (STYPE_INTEGER) 
      rdest%ivalue = int(dsource)
    case (STYPE_DOUBLE) 
      rdest%dvalue = dsource
    case (STYPE_INVALID) 
      ! Invalid and no variable associated? Will be a double.
      if (rdest%svarname .eq. "") then
        rdest%ctype = STYPE_DOUBLE
        rdest%dvalue = dsource
      else
        rdest%ctype = STYPE_INVALID
      end if
    case default
      ! Invalid value
      rdest%ctype = STYPE_INVALID
    end select
    
  end subroutine

  subroutine tpsym_setstring (rdest,ssource)
  type (t_symbolValue), intent(inout) :: rdest
  character(len=*), intent(in) :: ssource
  
    ! Type conversion if possible
    select case (rdest%ctype)
    case (STYPE_STRING) 
      rdest%svalue = ssource
      rdest%ilength = len_trim(ssource)
    case (STYPE_INVALID) 
      ! Invalid and no variable associated? Will be a double.
      if (rdest%svarname .eq. "") then
        rdest%ctype = STYPE_STRING
        rdest%svalue = ssource
        rdest%ilength = len_trim(ssource)
      else
        rdest%ctype = STYPE_INVALID
      end if
    case default
      ! Invalid value
      rdest%ctype = STYPE_INVALID
    end select
    
  end subroutine

  ! ***************************************************************************

  function tpsym_add (rsource1,rsource2) result (rdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  type (t_symbolValue) :: rdest
    
    rdest%ctype = STYPE_INVALID

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        rdest%ctype = STYPE_INTEGER
        rdest%ivalue = rsource1%ivalue + rsource2%ivalue
      
      case (STYPE_DOUBLE)
        rdest%ctype = STYPE_DOUBLE
        rdest%dvalue = real(rsource1%ivalue,DP) + rsource2%dvalue
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        rdest%ctype = STYPE_DOUBLE
        rdest%ivalue = rsource1%dvalue + real(rsource2%ivalue,DP)
      
      case (STYPE_DOUBLE)
        rdest%ctype = STYPE_DOUBLE
        rdest%dvalue = rsource1%dvalue + rsource2%dvalue
      end select
      
    case (STYPE_STRING)
      select case (rsource2%ctype)
      case (STYPE_STRING)
        ! Concatenate
        rdest%ctype = STYPE_STRING
        rdest%ilength = rsource1%ilength + rsource2%ilength
        rdest%svalue = rsource1%svalue(1:rsource1%ilength) // rsource2%svalue(1:rsource2%ilength)
      end select
    end select
    
  end function
  
  ! ***************************************************************************

  function tpsym_subtract (rsource1,rsource2) result (rdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  type (t_symbolValue) :: rdest
    
    rdest%ctype = STYPE_INVALID

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        rdest%ctype = STYPE_INTEGER
        rdest%ivalue = rsource1%ivalue - rsource2%ivalue
      
      case (STYPE_DOUBLE)
        rdest%ctype = STYPE_DOUBLE
        rdest%dvalue = real(rsource1%ivalue,DP) - rsource2%dvalue
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        rdest%ctype = STYPE_DOUBLE
        rdest%ivalue = rsource1%dvalue - real(rsource2%ivalue,DP)
      
      case (STYPE_DOUBLE)
        rdest%ctype = STYPE_DOUBLE
        rdest%dvalue = rsource1%dvalue - rsource2%dvalue
      end select
      
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_multiply (rsource1,rsource2) result (rdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  type (t_symbolValue) :: rdest
    
    rdest%ctype = STYPE_INVALID

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        rdest%ctype = STYPE_INTEGER
        rdest%ivalue = rsource1%ivalue * rsource2%ivalue
      
      case (STYPE_DOUBLE)
        rdest%ctype = STYPE_DOUBLE
        rdest%dvalue = real(rsource1%ivalue,DP) * rsource2%dvalue
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        rdest%ctype = STYPE_DOUBLE
        rdest%ivalue = rsource1%dvalue * real(rsource2%ivalue,DP)
      
      case (STYPE_DOUBLE)
        rdest%ctype = STYPE_DOUBLE
        rdest%dvalue = rsource1%dvalue * rsource2%dvalue
      end select
      
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_divide (rsource1,rsource2) result (rdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  type (t_symbolValue) :: rdest
    
    rdest%ctype = STYPE_INVALID

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        if (rsource2%ivalue .ne. 0) then
          rdest%ctype = STYPE_INTEGER
          rdest%ivalue = rsource1%ivalue / rsource2%ivalue
        end if
      
      case (STYPE_DOUBLE)
        if (rsource2%dvalue .ne. 0) then
          rdest%ctype = STYPE_DOUBLE
          rdest%dvalue = real(rsource1%ivalue,DP) / rsource2%dvalue
        end if
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        if (rsource2%ivalue .ne. 0) then
          rdest%ctype = STYPE_DOUBLE
          rdest%ivalue = rsource1%dvalue / real(rsource2%ivalue,DP)
        end if
      
      case (STYPE_DOUBLE)
        if (rsource2%dvalue .ne. 0) then
          rdest%ctype = STYPE_DOUBLE
          rdest%dvalue = rsource1%dvalue / rsource2%dvalue
        end if
      end select
      
    end select
    
  end function

! ***************************************************************************

  function tpsym_modulo (rsource1,rsource2) result (rdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  type (t_symbolValue) :: rdest
    
    rdest%ctype = STYPE_INVALID

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        if (rsource2%ivalue .ne. 0) then
          rdest%ctype = STYPE_INTEGER
          rdest%ivalue = mod(rsource1%ivalue, rsource2%ivalue)
        end if
      end select
      
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_eq (rsource1,rsource2) result (bdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  logical :: bdest
    
    ! Output is =0: false.
    ! May be set to =1: true.
    bdest = .false.

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%ivalue .eq. rsource2%ivalue)
      
      case (STYPE_DOUBLE)
        bdest =  (real(rsource1%ivalue,dp) .eq. rsource2%dvalue)
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%dvalue .eq. real(rsource2%ivalue,dp))
      
      case (STYPE_DOUBLE)
        bdest =  (rsource1%dvalue .eq. rsource2%dvalue)
      end select
      
    case (STYPE_STRING)
      select case (rsource2%ctype)
      
      case (STYPE_STRING)
        bdest =  (rsource1%svalue .eq. rsource2%svalue)
      end select
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_neq (rsource1,rsource2) result (bdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  logical :: bdest
    
    ! Output is =0: false.
    ! May be set to =1: true.
    bdest = .false.

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%ivalue .ne. rsource2%ivalue)
      
      case (STYPE_DOUBLE)
        bdest =  (real(rsource1%ivalue,dp) .ne. rsource2%dvalue)
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%dvalue .ne. real(rsource2%ivalue,dp))
      
      case (STYPE_DOUBLE)
        bdest =  (rsource1%dvalue .ne. rsource2%dvalue)
      end select
      
    case (STYPE_STRING)
      select case (rsource2%ctype)
      
      case (STYPE_STRING)
        bdest =  (rsource1%svalue .ne. rsource2%svalue)
      end select
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_less (rsource1,rsource2) result (bdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  logical :: bdest
    
    ! Output is =0: false.
    ! May be set to =1: true.
    bdest = .false.

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%ivalue .lt. rsource2%ivalue)
      
      case (STYPE_DOUBLE)
        bdest =  (real(rsource1%ivalue,dp) .lt. rsource2%dvalue)
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%dvalue .lt. real(rsource2%ivalue,dp))
      
      case (STYPE_DOUBLE)
        bdest =  (rsource1%dvalue .lt. rsource2%dvalue)
      end select
      
    case (STYPE_STRING)
      select case (rsource2%ctype)
      
      case (STYPE_STRING)
        bdest =  (rsource1%svalue .lt. rsource2%svalue)
      end select
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_lesseq (rsource1,rsource2) result (bdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  logical :: bdest
    
    ! Output is =0: false.
    ! May be set to =1: true.
    bdest = .false.

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%ivalue .le. rsource2%ivalue)
      
      case (STYPE_DOUBLE)
        bdest =  (real(rsource1%ivalue,dp) .le. rsource2%dvalue)
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%dvalue .le. real(rsource2%ivalue,dp))
      
      case (STYPE_DOUBLE)
        bdest =  (rsource1%dvalue .le. rsource2%dvalue)
      end select
      
    case (STYPE_STRING)
      select case (rsource2%ctype)
      
      case (STYPE_STRING)
        bdest =  (rsource1%svalue .le. rsource2%svalue)
      end select
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_greater (rsource1,rsource2) result (bdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  logical :: bdest
    
    ! Output is =0: false.
    ! May be set to =1: true.
    bdest = .false.

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%ivalue .gt. rsource2%ivalue)
      
      case (STYPE_DOUBLE)
        bdest =  (real(rsource1%ivalue,dp) .gt. rsource2%dvalue)
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%dvalue .gt. real(rsource2%ivalue,dp))
      
      case (STYPE_DOUBLE)
        bdest =  (rsource1%dvalue .gt. rsource2%dvalue)
      end select
      
    case (STYPE_STRING)
      select case (rsource2%ctype)
      
      case (STYPE_STRING)
        bdest =  (rsource1%svalue .gt. rsource2%svalue)
      end select
    end select
    
  end function

  ! ***************************************************************************

  function tpsym_greatereq (rsource1,rsource2) result (bdest)
  type (t_symbolValue), intent(in) :: rsource1
  type (t_symbolValue), intent(in) :: rsource2
  logical :: bdest
    
    ! Output is =0: false.
    ! May be set to =1: true.
    bdest = .false.

    ! The type decides. Maybe there must be a type conversion.
    select case (rsource1%ctype)
    case (STYPE_INTEGER)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%ivalue .ge. rsource2%ivalue)
      
      case (STYPE_DOUBLE)
        bdest =  (real(rsource1%ivalue,dp) .ge. rsource2%dvalue)
      end select

    case (STYPE_DOUBLE)
      select case (rsource2%ctype)
      
      case (STYPE_INTEGER)
        bdest =  (rsource1%dvalue .ge. real(rsource2%ivalue,dp))
      
      case (STYPE_DOUBLE)
        bdest =  (rsource1%dvalue .ge. rsource2%dvalue)
      end select
      
    case (STYPE_STRING)
      select case (rsource2%ctype)
      
      case (STYPE_STRING)
        bdest =  (rsource1%svalue .ge. rsource2%svalue)
      end select
    end select
    
  end function

end module
