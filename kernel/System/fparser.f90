!##############################################################################
!# ****************************************************************************
!# <name> fparser </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This public domain function parser module is intended for applications 
!# where a set of mathematical expressions is specified at runtime and is 
!# then evaluated for a large number of variable values. This is done by
!# compiling the set of function strings into byte code, which is interpreted
!# very efficiently for the various variable values. The evaluation is 
!# straightforward and no recursions are done (uses stack arithmetic).
!#
!# ------------------------------------------------------------------------
!# Copyright notice
!# ------------------------------------------------------------------------
!# (a) This module is based on the "Fortran 90 function parser V1.1" 
!#     written by Roland Schmehl < Roland.Schmehl@mach.uni-karlsruhe.de >
!#     The original Fortran90 source code is available from:
!#     http://itsextern.its.uni-karlsruhe.de/~schmehl/functionparser.html
!#
!#     However the "Fortran 90 function parser V1.1" only recognizes the 
!#     (single argument) Fortran 90 intrinsic functions abs, exp, log10,
!#     log, sqrt, sinh, cosh, tanh, sin, cos, tan, asin, acos, atan.
!#
!#     The function parser concept is based on a C++ class library written
!#     by Warp < warp@iki.fi > available from:
!#     http://warp.povusers.org/FunctionParser/
!#
!# (b) This FParser module is an extension of the "Fortran 90 function parser
!#     V1.1" which implements most of the features available in the "function
!#     parser library for C++ V2.8" written by Warp. The optimizer included
!#     in the C++ library and the recursive evaluation of functions by means
!#     of eval(...) is not implemented in this version.
!#
!# ------------------------------------------------------------------------
!# Basic usage
!# ------------------------------------------------------------------------
!#
!# Step 0 - Module Import
!# ----------------------
!# In all program units where you want to use the function parser procedures 
!# and variables you must import the module by:
!#
!# USE fparser
!#
!# This command imports only 6 public names: fparser_create, fparser_release, 
!# fparser_parseFunction, fparser_evalFunction, fparser_ErrorMsg and EvalErrType
!# which are explained in the following. The remainder of the 
!# module is hidden to the calling program.
!#
!# Step 1 - Initialization
!# -----------------------
!# The parser module has to be initialized for the simultaneous evaluation of 
!# n functions by calling the module subroutine initp one time in your Fortran 
!# code:
!# 
!# CALL fparser_create (Parser, n)
!# 
!# This allocates i=1,...,n internal data structures used by the byte-compiler 
!# and subsequently by the bytecode-interpreter in the bytecode object Comp.
!#
!# Step 2 - Function parsing
!# -------------------------
!# The i-th function string FuncStr is parsed (checked and compiled) into the 
!# i-th bytecode by calling the module subroutine parsef:
!#
!# CALL fparser_parseFunction (Parser, i, FuncStr, Var)
!#
!# The variable names as they appear in the string FuncStr have to be passed 
!# in the one-dimensional string array Var (zero size of Var is acceptable). 
!# The number of variables is implicitly passed by the dimension of this array. 
!# For some notes on the syntax of the function string see below.
!#
!# Step 3 - Function evaluation
!# ----------------------------
!# The i-th function value is evaluated for a specific set of variable values 
!# by calling the module function evalf:
!#
!# a = fparser_evalFunction (Parser, i, Val)
!#
!# The variable values are passed in the one-dimensional array Val which must 
!# have the same dimension as array Var. 
!#
!# ------------------------------------------------------------------------
!# Error handling
!# ------------------------------------------------------------------------
!# 
!# An error in the function parsing step leads to a detailed error message 
!# (Type and position of error) and program termination.
!#
!# An error during function evaluation returns a function value of 0.0 and
!# sets the error flag EvalErrType (part of the t_fparser derived type) to 
!# a value > 0 (EvalErrType=0 indicates no error). An error message from the 
!# bytecode-interpreter can be obtained by calling the character function 
!# fparser_ErrorMsg (Parser) with the parser object as an argument.
!#
!# ------------------------------------------------------------------------
!# Function string syntax
!# ------------------------------------------------------------------------
!#
!# Although they have to be passed as array elements of the same declared 
!# length (Fortran 90 restriction), the variable names can be of arbitrary 
!# actual length for the parser. Parsing for variables is case sensitive. 
!#
!# The syntax of the function string is similar to the Fortran convention. 
!# Mathematical Operators recognized are +, -, *, /, %, ** or alternatively
!# ^, whereas symbols for brackets must be (), [] or {}. Note that the 
!# parser does not check if, e.g. ( is closed by ) or ]. At the moment,
!# different brackets may be used only to improve readability of the function
!# string.
!#
!# Operations are evaluated in the correct order:
!#
!#  ()             expressions in brackets first
!#  -A             unary minus (or plus)
!#  A**B A^B       exponentiation (A raised to the power B)
!#  A*B  A/B  A%B  multiplication, division and modulo
!#  A+B  A-B       addition and subtraction
!#  A=B  A!=B  A < B  A <= B  A > B  A >= B
!#                 comparison between A and B (result is either 0 or 1)
!#  A&B            result is 1 if int(A) and int(B) differ from 0, else 0.
!#  A|B            result is 1 if int(A) or int(B) differ from 0, else 0.
!#
!# The function string can contain integer or real constants. To be recognized
!# as explicit constants these must conform to the format
!#
!# [+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
!#
!# where nnn means any number of digits. The mantissa must contain at least
!# one digit before or following an optional decimal point. Valid exponent 
!# identifiers are 'e', 'E', 'd' or 'D'. If they appear they must be followed 
!# by a valid exponent!
!#
!# Note that the function parser is case insensitive.
!# The following mathematical functions are supported
!#
!# abs(A)    : Absolute value of A. If A is negative, returns -A otherwise
!#             returns A.
!# acos(A)   : Arc-cosine of A. Returns the angle, measured in radians,
!#             whose cosine is A.
!# acosh(A)  : Same as acos() but for hyperbolic cosine.
!# aint(A)   : Truncate A to a whole number
!# anint(A)  : Rounds A to the closest integer. 0.5 is rounded to 1.
!# asin(A)   : Arc-sine of A. Returns the angle, measured in radians, whose
!#             sine is A.
!# asinh(A)  : Same as asin() but for hyperbolic sine.
!# atan(A)   : Arc-tangent of (A). Returns the angle, measured in radians,
!#             whose tangent is (A).
!# atan2(A,B): Arc-tangent of A/B. The two main differences to atan() is
!#             that it will return the right angle depending on the signs of
!#             A and B (atan() can only return values betwen -pi/2 and pi/2),
!#             and that the return value of pi/2 and -pi/2 are possible.
!# atanh(A)  : Same as atan() but for hyperbolic tangent.
!# ceil(A)   : Ceiling of A. Returns the smallest integer greater than A.
!#             Rounds up to the next higher integer.
!# cos(A)    : Cosine of A. Returns the cosine of the angle A, where A is
!#             measured in radians.
!# cosh(A)   : Same as cos() but for hyperbolic cosine.
!# cot(A)    : Cotangent of A (equivalent to 1/tan(A)).
!# csc(A)    : Cosecant of A (equivalent to 1/sin(A)).
!# exp(A)    : Exponential of A. Returns the value of e raised to the power
!#             A where e is the base of the natural logarithm, i.e. the
!#             non-repeating value approximately equal to 2.71828182846.
!# floor(A)  : Floor of A. Returns the largest integer less than A. Rounds
!#             down to the next lower integer.
!# if(A,B,C) : If int(A) differs from 0, the return value of this function is B,
!#             else C. Only the parameter which needs to be evaluated is
!#             evaluated, the other parameter is skipped.
!# log(A)    : Natural (base e) logarithm of A.
!# log10(A)  : Base 10 logarithm of A.
!# max(A,B)  : If A > B, the result is A, else B.
!# min(A,B)  : If A < B, the result is A, else B.
!# sec(A)    : Secant of A (equivalent to 1/cos(A)).
!# sin(A)    : Sine of A. Returns the sine of the angle A, where A is
!#             measured in radians.
!# sinh(A)   : Same as sin() but for hyperbolic sine.
!# sqrt(A)   : Square root of A. Returns the value whose square is A.
!# tan(A)    : Tangent of A. Returns the tangent of the angle A, where A
!#             is measured in radians.
!# tanh(A)   : Same as tan() but for hyperbolic tangent.
!#
!# The parser also supports a number of standard constants in the function
!# string. All constants start with an underscore '_'. The following constants
!# are defined by default:
!#
!# _PI       : Gives the number $pi$.
!# _EXP      : Gives the number $e$
!# _INFTY    : Gives the maximum possible number in double precision, defined
!#             in fsystem by SYS_INFINITY.
!#
!# In addition, the use can define his own global constant which are available
!# throughout the complete function parser as it is the case for the standard
!# constant defined above.
!#
!# The parser also supports user-defined expressions which are globally
!# available throughout the complete function parser. All expressions must
!# start with '@' to indicate that the following expression should be
!# looked-up from the list of predefined expressions.
!#
!# The following routines can be found in this module:
!#
!# 1.) fparser_init
!#     -> Initialize the sub-system for function parsers
!#
!# 2.) fparser_done
!#     -> Release the sub-system for function parsers
!#
!# 3.) fparser_defineConstant
!#     -> Define special constants which are available for all function parsers
!#
!# 4.) fparser_defineExpression
!#     -> Define special expressions which are available for all function parsers
!#
!# 5.) fparser_create
!#     -> Create function parser
!#
!# 6.) fparser_release
!#     -> Release function parser
!#
!# 7.) fparser_parseFunction
!#     -> Parse function string and compile it into bytecode
!#
!# 8.) fparser_evalFunction = fparser_evalFunctionScalar /
!#                            fparser_evalFunctionArray
!#     -> Evaluate bytecode
!#
!# 9.) fparser_ErrorMsg
!#     -> Get error message from function parser
!#
!# 10.) fparser_PrintByteCode
!#      -> Print the bytecode stack (very technical!)
!#
!# 11.) fparser_parseFileForKeyword
!#      -> Parse input file for keyworkd
!#
!# The following internal routines can be found in this module:
!#
!# 1.) CheckSyntax
!#     -> Check syntax of function string before compiling bytecode
!#
!# 2.) isOperator
!#     -> Return size of operator and 0 otherwise
!#
!# 3.) MathFunctionIndex
!#     -> Return index of mathematical function and 0 otherwise
!#
!# 4.) MathFunctionParameters
!#     -> Return number of required function parameters
!#
!# 5.) ConstantIndex
!#     -> Return index of predefined constant and 0 otherwise
!#
!# 6.) ExpressionIndex
!#     -> Return index of predefined expression and 0 otherwise
!#
!# 7.) VariableIndex
!#     -> Return index of variable
!#
!# 8.) RemoveSpaces
!#     -> Remove spaces from string
!#
!# 9.) Replace
!#     -> Replace all appearances of one character set by another
!#        character set in a given string
!#
!# 10.) Compile
!#      -> Compile function string into bytecode
!#
!# 11.) incStackPtr
!#     -> Increase stack pointer
!#
!# 12.) AddCompiledByte
!#     -> Add compiled byte to bytecode stack
!#
!# 13.) RemoveCompiledByte
!#     -> Remove last compiled byte from bytecode stack
!#
!# 14.) AddImmediate
!#     -> Add immediate to immediate stack
!#
!# 15.) AddFunctionOpcode
!#     -> Add function opcode to bytecode stack
!#
!# 16.) RealNum
!#     -> Get real number from string
!#
!# 17.) FunctionSize
!#     -> Get the total size of the function
!#
!# 18.) CompileExpression
!#      -> Compile ','
!#
!# 19.) CompileOr
!#      -> Compile '|'
!#
!# 20.) CompileAnd
!#      -> Compile '&'
!#
!# 21.) CompileComparison
!#      -> Compile '=', '<', and '>'
!#
!# 22.) CompileAddition
!#      -> Compile '+' and '-'
!#
!# 23.) CompileMult
!#      -> Compile '*', '/', and '%'
!#
!# 24.) CompileUnaryMinus
!#      -> Compile unary '-'
!#
!# 25.) CompilePow
!#      -> Compile '^'
!#
!# 26.) CompileElement
!#      -> Compile mathematical function, variable, constant and number
!#
!# 27.) CompileFunctionParameters
!#      -> Compile function parameters
!#
!# 28.) CompileIf
!#      -> Compile if-then-else
!#
!# </purpose>
!##############################################################################

module fparser
  
  use fsystem
  use genoutput
  use io
  use storage
  use stack
  
  implicit none
  
  private
  public :: t_fparser
  public :: fparser_init
  public :: fparser_done
  public :: fparser_parseFileForKeyword
  public :: fparser_defineConstant
  public :: fparser_defineExpression
  public :: fparser_create
  public :: fparser_release
  public :: fparser_parseFunction
  public :: fparser_evalFunction
  public :: fparser_ErrorMsg
  public :: fparser_PrintByteCode
  
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  
  interface fparser_evalFunction
    module procedure fparser_evalFunctionScalar
    module procedure fparser_evalFunctionBlock
  end interface

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<constants>

!<constantblock description="Global constants for parser">
  ! Maximum stack size
  ! The stack is used for all instances of function parsers. During initialization
  ! a smaller size can be specified. This is only the maximum size of the stack.
  integer, parameter :: FPAR_MAXSTACKSIZE   = 1024*1024 ! this is 8 MByte

  ! Length of string
  integer, parameter :: FPAR_STRLEN         = 2048

  ! Maximum number of predefined/user-defined constants
  integer, parameter :: FPAR_MAXCONSTS      = 128

  ! Maximum number of predefined/user-defined expressions
  integer, parameter :: FPAR_MAXEXPRESSIONS = 128

  ! Length of constant name
  integer, parameter :: FPAR_CONSTLEN       = 12

  ! Length of expression name
  integer, parameter :: FPAR_EXPRLEN        = 12
!</constantblock>


!<constantblock description="types for parser expressions">
  
  ! Constant
  integer, parameter, public :: FPAR_CONSTANT   = 1

  ! Expression
  integer, parameter, public :: FPAR_EXPRESSION = 2
!</constantblock>


!<constantblock description="data type for parser bytecode">
  ! Data type of bytecode
  integer, parameter :: is = selected_int_kind(1)
!</constantblock>


!<constantblock description="keywords for parser">
  integer(is), parameter :: cImmed       =  1, &
                            cJump        =  2, &
                            cNeg         =  3, &
                            cDeg         =  4, &
                            cRad         =  5, &
                            cAdd         =  6, & ! <-- first dyadic operator: A .OP. B
                            cSub         =  7, &
                            cMul         =  8, &
                            cDiv         =  9, &
                            cMod         = 10
  integer(is), parameter :: cPow         = 11, &
                            cNEqual      = 12, & ! NOTE: != must be prior to =
                            cEqual       = 13, &
                            cLessOrEq    = 14, & ! NOTE: <= must be prior to < 
                            cLess        = 15, &
                            cGreaterOrEq = 16, & ! NOTE: >= must be prior to > 
                            cGreater     = 17, &
                            cNot         = 18, &
                            cAnd         = 19, &
                            cOr          = 20    ! --> last dyadic operator: A.OP. B
  integer(is), parameter :: cIf          = 21, & ! <-- if-then-else
                            cMin         = 22, & ! <-- first dyadic operator: .OP.(A,B)
                            cMax         = 23, &
                            cAtan2       = 24, & ! --> last dyadic operator: .OP.(A,B)
                            cAbs         = 25, & ! <-- monadic operator: .OP.(A)
                            cAnint       = 26, &
                            cAint        = 27, &
                            cExp         = 28, &
                            cLog10       = 29, &
                            cLog         = 30
  integer(is), parameter :: cSqrt        = 31, &
                            cSinh        = 32, &
                            cCosh        = 33, &
                            cTanh        = 34, &
                            cSin         = 35, &
                            cCos         = 36, &
                            cTan         = 37, & 
                            cCot         = 38, &
                            cAsinh       = 39, &
                            cAsin        = 40
  integer(is), parameter :: cAcosh       = 41, &
                            cAcos        = 42, &
                            cAtanh       = 43, &
                            cAtan        = 44, &
                            cCeil        = 45, &
                            cFloor       = 46, &
                            cCsc         = 47, &
                            cSec         = 48, & ! --> last monadic operator: .OP.(A)
                            VarBegin     = 49
!</constantblock>


!<constantblock description="symbols for parser operands">
  character (LEN=2), dimension(cAdd:cOr), parameter :: Ops = (/ '+ ', &
                                                                '- ', &
                                                                '* ', &
                                                                '/ ', &
                                                                '% ', &
                                                                '^ ', &
                                                                '!=', &
                                                                '= ', &
                                                                '<=', &
                                                                '< ', &
                                                                '>=', &
                                                                '> ', &
                                                                '! ', &
                                                                '& ', &
                                                                '| ' /)
!</constantblock>


!<constantblock description="function names for parser">
  character (LEN=5), dimension(cIf:cSec), parameter :: Funcs = (/ 'if   ', &
                                                                  'min  ', &
                                                                  'max  ', &
                                                                  'atan2', &
                                                                  'abs  ', &
                                                                  'anint', &
                                                                  'aint ', &
                                                                  'exp  ', &
                                                                  'log10', &
                                                                  'log  ', &
                                                                  'sqrt ', &
                                                                  'sinh ', &
                                                                  'cosh ', &
                                                                  'tanh ', &
                                                                  'sin  ', &
                                                                  'cos  ', &
                                                                  'tan  ', &
                                                                  'cot  ', &
                                                                  'asinh', &
                                                                  'asin ', &
                                                                  'acosh', &
                                                                  'acos ', &
                                                                  'atanh', &
                                                                  'atan ', &
                                                                  'ceil ', &
                                                                  'floor', &
                                                                  'csc  ', &
                                                                  'sec  ' /)
!</constantblock>


!<constantblock description="predefined constant names for parser; an underscore '_' is automatically added">
  character(LEN=FPAR_CONSTLEN), dimension(3) :: PredefinedConsts = (/ 'pi        ', &
                                                                      'exp       ', &
                                                                      'infty     ' /)
!</constantblock>


!<constantblock description="predefined constant values for parser">
  real(DP), dimension(3), parameter :: PredefinedConstvals = (/&
      3.141592653589793115997963468544185161590576171875_DP, &
      2.718281828459045090795598298427648842334747314453125_DP, &
      SYS_INFINITY/)
!</constantblock>

!</constants>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<globals>

  ! Global handle to stack memory
  integer, save                                                    :: h_Stack = ST_NOHANDLE
  
  ! Global pointer to stack memory
  real(DP), dimension(:), pointer, save                            :: p_Stack => null()

  ! Global number of predefined/user-defined constants
  integer, save                                                    :: nConsts = 0 

  ! Global constant names for parser
  character(LEN=FPAR_CONSTLEN), dimension(FPAR_MAXCONSTS), save    :: Consts  = '     '

  ! Global constant values for parser
  real(DP), dimension(FPAR_MAXCONSTS), save                        :: ConstVals = 0

  ! Global number of predefined/user-defined expressions
  integer, save                                                    :: nExpressions = 0

  ! Global expression name for parser
  character(LEN=FPAR_EXPRLEN), dimension(FPAR_MAXCONSTS), save     :: Expressions = ''

  ! Global expression string for parser
  character(LEN=FPAR_STRLEN), dimension(FPAR_MAXEXPRESSIONS), save :: ExpressionStrings
!</globals>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<types>

!<typeblock>
  
  ! Type block for storing all information of the function parser
  type t_fparser
    private
    
    ! Array of function parser components. Each component is used to
    ! handle one function string at a time
    type(t_fparserComponent), dimension(:), pointer :: Comp => null()
    
    ! Number of parser components
    integer :: nComp                                        = 0
  end type t_fparser

!</typeblock>

!<typeblock>

  ! Type block for storing the bytecode of the function parser for one component
  type t_fparserComponent
    private

    ! Size of bytecode
    integer :: ByteCodeSize                                 = 0

    ! Size of immediates
    integer :: ImmedSize                                    = 0

    ! Stack size
    integer :: StackSize                                    = 0

    ! Stack pointer
    integer :: StackPtr                                     = 0

    ! Use degree conversion DEG <-> RAD for some functions
    logical :: useDegreeConversion                          = .false.

    ! Is vectorizable
    logical :: isVectorizable                               = .true.
    
    ! Bytecode
    integer(is), dimension(:), pointer :: ByteCode          => null()
    
    ! Immediates
    real(DP), dimension(:), pointer :: Immed                => null()
  end type t_fparserComponent
!</typeblock>

!</types>

contains

  ! *****************************************************************************
  
!<subroutine>

  subroutine fparser_init (istacksize)

!<description>
    ! Initialize function parser
!</description>

!<input>
    ! OPTIONAL: initial size of the stack memory
    integer, intent(IN), optional :: istacksize
!</input>
!</subroutine>

    integer :: isize,isize1,iconst

    if (present(istacksize)) then
      isize=min(istacksize,FPAR_MAXSTACKSIZE)
    else
      isize=FPAR_MAXSTACKSIZE
    end if
    
    if (h_Stack .eq. ST_NOHANDLE) then
      
      ! Allocate memory for global stack and set pointer
      call storage_new('fparser_init','p_Stack',isize,ST_DOUBLE,h_Stack,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double(h_Stack,p_Stack)

      ! Initialize predefined constants
      do iconst=lbound(PredefinedConsts,1),ubound(PredefinedConsts,1)
        call fparser_defineConstant(PredefinedConsts(iconst),PredefinedConstvals(iconst))
      end do

    else
      
      ! Realloc memory for global stack if required
      call storage_getsize(h_Stack,iSize1)
      if (isize1 .lt. isize) then
        call storage_realloc('fparser_init',isize,h_Stack,ST_NEWBLOCK_NOINIT,.true.)
        call storage_getbase_double(h_Stack,p_Stack)
      end if

    end if   
    
  end subroutine fparser_init

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_done

!<description>
    ! Release function parser
!</description>
!</subroutine>

    ! Reset constants
    nConsts   = 0
    Consts    = '     '
    ConstVals = 0._DP
    
    ! Reset expressions
    nExpressions      = 0
    Expressions       = ''
    ExpressionStrings = ''

    ! Free memory
    if (h_Stack  /= ST_NOHANDLE) call storage_free(h_Stack)
    nullify(p_Stack)
  end subroutine fparser_done

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_parseFileForKeyword(sfilename, ckeyword, itype)

!<description>
    ! Parse the file for the given keyword and make it a constant or a
    ! predefined expression depending on the variable itype
!</description>

!<input>
    ! name of parameter file to be parsed
    character(LEN=*), intent(IN) :: sfilename

    ! name of keyword to parser for
    character(LEN=*), intent(IN) :: ckeyword

    ! type of keyword: FPAR_CONSTANT, FPAR_EXPRESSION
    integer, intent(IN) :: itype
!</input>
!</subroutine>

    ! local variables
    character(SYS_STRLEN)  :: skeyword,sname
    character(FPAR_STRLEN) :: sdata,svalue
    real(DP) :: dvalue
    integer :: iunit,ipos,jpos,ios,idatalen

    ! Try to open the file
    call io_openFileForReading(sfilename, iunit, .true.)

    ! Oops...
    if (iunit .eq. -1) then
      call output_line('Unable to open input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'fparser_parseFileForKeyword')
      call sys_halt()
    end if

    ! Read through the complete input file and look for global definitions
    ! of constants and fixed expressions
    ios = 0
    do while(ios .eq. 0)

      ! Read next line in file
      call io_readlinefromfile(iunit, sdata, idatalen, ios)
      
      ! Check for keyword defconst or defexp
      ipos = scan(sdata(1:idatalen), ":")
      if (ipos .eq. 0) cycle
      
      call sys_tolower(sdata(1:max(1,ipos-1)), skeyword)
      if (trim(adjustl(skeyword)) .eq. ckeyword) then
        
        ! Split the line into name and value
        jpos  = scan(sdata(1:idatalen), "=" , .true.)
        sname  = trim(adjustl(sdata(ipos+1:jpos-1)))
        svalue = trim(adjustl(sdata(jpos+1:)))
                
        ! We found a keyword that will be applied to the parser
        select case(itype)
        case(FPAR_CONSTANT)
          read(svalue,*) dvalue
          call fparser_defineConstant(sname, dvalue)

        case(FPAR_EXPRESSION)
          call fparser_defineExpression(sname, svalue)

        case DEFAULT
          call output_line('Invalid type of expression!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'fparser_parseFileForKeyword')
          call sys_halt()
        end select

      end if
    end do

    ! Close file
    close (iunit)

  end subroutine fparser_parseFileForKeyword
  
  ! *****************************************************************************

!<subroutine>

  subroutine fparser_defineConstant(constant,constantValue)

!<description>
    ! Define a new constant for the function parser.
    ! This subroutine checks if the given constant is already defined.
!</description>

!<input>
    ! Name of the constant
    character(LEN=FPAR_CONSTLEN), intent(IN) :: constant

    ! Value of the constant
    real(DP), intent(IN) :: constantValue
!</input>
!</subroutine>

    ! local variables
    character(LEN=len(constant)) :: const
    integer :: iConst

    ! Check if there is enough space
    if (nConsts < FPAR_MAXCONSTS) then
      nConsts = nConsts+1
    else
      call output_line('No space left for definition of constant!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'fparser_defineConstant')
      call sys_halt()
    end if

    ! Prepare constant
    call sys_tolower(constant,const)

    ! Check if constant is already defined
    do iConst=1,nConsts-1
      if (Consts(iConst) == const) then
        ! If it is already defined, then it must not have a different value
        if(ConstVals(iConst) /= constantValue) then
          call output_line('Constant is already defined with different value!',&
                            OU_CLASS_WARNING,OU_MODE_STD,'fparser_defineConstant')
          call sys_halt()
        else
          nConsts = nConsts-1
          return
        end if
      end if
    end do

    ! Apply constant value and constant name
    Consts(nConsts)    = const
    ConstVals(nConsts) = constantValue
  end subroutine fparser_defineConstant

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_defineExpression(expression,expressionString)

!<description>
    ! Define a new expression for the function parser.
    ! This subroutine checks if the given expression is already defined.
!</description>

!<input>
    ! Name of the expression
    character(LEN=FPAR_EXPRLEN), intent(IN) :: expression

    ! String of the expression
    character(LEN=*), intent(IN) :: expressionString
!</input>
!</subroutine>

    ! local variables
    character(LEN=len(expression))       :: expr
    character(LEN=len(expressionString)) :: str
    integer :: iExpression

    ! Check if there is enough space
    if (nExpressions < FPAR_MAXEXPRESSIONS) then
      nExpressions = nExpressions+1
    else
      call output_line('No space left for definition of expression!',&
                        OU_CLASS_WARNING,OU_MODE_STD,'fparser_defineExpression')
      call sys_halt()
    end if

    ! Prepare expression string
    call sys_tolower(expression,expr)
    call sys_tolower(expressionString,str)

    ! Replace human readable function names by 1-Char. format
    call Replace ('**','^ ',str)
    call Replace ('[','(',  str)
    call Replace (']',')',  str)
    call Replace ('{','(',  str)
    call Replace ('}',')',  str)
    
    ! Condense function string
    call RemoveSpaces (str)

    ! Check if expressions is already defined
    do iExpression=1,nExpressions-1
      if (Expressions(iExpression) == expr) then
        ! If it is already defined, then it must not have a different value
        if(ExpressionStrings(iExpression) /= str) then
          call output_line('Expression is already defined with different string!',&
                            OU_CLASS_WARNING,OU_MODE_STD,'fparser_defineExpression')
          call sys_halt()
        else
          nExpressions = nExpressions-1
          return
        end if
      end if
    end do

    ! Apply expressions string and expression name
    Expressions(nExpressions)       = expr
    ExpressionStrings(nExpressions) = str
  end subroutine fparser_defineExpression

  ! *****************************************************************************

!<subroutine>
  
  subroutine fparser_create (rparser,nComp)

!<description>
    ! Initialize function parser for nComp functions
!</description>

!<input>
    ! Number of functions
    integer, intent(IN) :: nComp
!</input>

!<inputoutput>
    ! Function parser object
    type(t_fparser), intent(INOUT) :: rparser
!</inputoutput>
!</subroutine>
    
    ! Set number of components
    rparser%nComp=nComp
    
    if (associated(rparser%Comp)) deallocate(rparser%Comp)
    allocate (rparser%Comp(nComp))
  end subroutine fparser_create

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_release (rparser)

!<description>
    ! Release function parser and all of its coponents
!</description>

!<inputoutput>
    ! Function parser
    type(t_fparser), intent(INOUT) :: rparser
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: iComp
    
    ! Check that pointer is associated and return otherwise
    if (.not.associated(rparser%Comp)) return

    ! Loop over all components
    do iComp=1,rparser%nComp
      if (associated(rparser%Comp(iComp)%ByteCode)) deallocate(rparser%Comp(iComp)%ByteCode)
      if (associated(rparser%Comp(iComp)%Immed))    deallocate(rparser%Comp(iComp)%Immed)
    end do
    deallocate(rparser%Comp)
    rparser%nComp=0
  end subroutine fparser_release

  ! *****************************************************************************

!<subroutine>
  
  subroutine fparser_parseFunction (rparser, iComp, FuncStr, Var, useDegrees)

!<description>
    ! Parse ith function string FuncStr and compile it into bytecode
!</description>
    
!<input>
    ! Function identifier
    integer, intent(IN) :: iComp

    ! Function string
    character (LEN=*), intent(IN) :: FuncStr

    ! Array with variable names
    character (LEN=*), dimension(:), intent(IN) :: Var

    ! OPTIONAL
    logical, intent(IN), optional :: useDegrees
!</input>

!<inputoutput>
    ! Function parser
    type (t_fparser), intent(INOUT) :: rparser
!</inputoutput>
!</subroutine>

    ! local variables
    character (LEN=len(FuncStr)) :: Func
    
    if (iComp < 1 .or. iComp > rparser%nComp) then
      write(*,*) '*** Parser error: Function number ',iComp,' out of range'
      call sys_halt()
    end if
    
    ! Local copy of function string
    Func = FuncStr

    ! Replace human readable function names by 1-Char. format
    call Replace ('**','^ ',Func)
    call Replace ('[','(',Func)
    call Replace (']',')',Func)
    call Replace ('{','(',Func)
    call Replace ('}',')',Func)

    ! Condense function string
    call RemoveSpaces (Func)

    ! Check for valid syntax; this prevents the bytecode compiler
    ! from running into endless loops or other problems
    call CheckSyntax (Func,Var)
    
    ! If syntax is correct, then compile into bytecode
    if (present(useDegrees))&
        rparser%Comp(iComp)%useDegreeConversion = useDegrees

    call Compile (rparser%Comp(iComp),Func,Var)
  end subroutine fparser_parseFunction

  ! *****************************************************************************

!<subroutine>
  
  subroutine fparser_evalFunctionScalar (rparser, iComp, Val, Res)

!<description>
    ! Evaluate bytecode of component iComp for the values passed in array 
    ! Val(:). Note that this function is a wrapper for the working routine
    ! evalFunctionScalar. It is used to adjust the dimensions of the global
    ! stack memory if required.
!</description>

!<input>
    ! Function parser
    type (t_fparser),  intent(IN) :: rparser

    ! Function identifier
    integer, intent(IN) :: iComp

    ! Variable values
    real(DP), dimension(:), intent(IN) :: Val
!</input>

!<output>
    ! Evaluated function
    real(DP), intent(OUT)  :: Res
!</output>
!</subroutine>

    ! local variables
    integer :: EvalErrType

    if (h_Stack .eq. ST_NOHANDLE) then
      call output_line('Parser not initialised!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'fparser_evalFunctionScalar')
      call sys_halt()
    end if
    
    ! Check if memory of global stack is sufficient
    if (size(p_Stack) < rparser%Comp(iComp)%StackSize+1) then
      if (rparser%Comp(iComp)%StackSize+1 < FPAR_MAXSTACKSIZE) then
        call storage_realloc('fparser_evalFunctionScalar',&
            rparser%Comp(iComp)%StackSize+1,h_Stack,ST_NEWBLOCK_NOINIT,.false.)
        call storage_getbase_double(h_Stack,p_Stack)
      else
        call output_line('Stack size exceeds memory!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'fparser_evalFunctionScalar')
        call sys_halt()
      end if
    end if
    
    ! Invoke working routine
    call evalFunctionScalar(p_Stack,rparser%Comp(iComp),Val,EvalErrType,Res)

    ! Check if evaluation was successful
    if (EvalErrType .ne. 0) then
      call output_line('An error occured during function evaluation!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'fparser_evalFunctionScalar')
      call sys_halt()
    end if
  end subroutine fparser_evalFunctionScalar

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_evalFunctionBlock (rparser, iComp, iDim, ValBlock, Res, ValScalar)

!<description>
    ! Evaluate bytecode of component iComp for the array of values passed in 
    ! array ValBlock(:,:). Note that this function is a wrapper for the working 
    ! routine evalFunctionBlock. It is used to adjust the dimensions of the 
    ! global stack memory if required. In some situations, there a variables,
    ! such as nodal coordinates, which are different for each component of 
    ! the resulting vector and those which are the same, e.g., time variable.
    ! The latter ones can be passed to the ValScalar argument which is used
    ! uniformly for each component of Res.
    !
    ! WARNING: The ordering of the variables must be identical to that given
    ! during the byte-code compilation. Care must be taken by the user since
    ! this cannot be checked be the function parser. Hence, if both ValBlock
    ! and ValScalar should be used, then the first variables must stored as
    ! blocks whereas the last variables can be scalar. This sound slightly
    ! complicated but here is an example:
    !
    ! Suppuse you want to evaluate a function f=f(x,y,t). You know that x,y
    ! corresponds to the coordinate vector and t denotes the time. Then
    ! you should order your variables according to [x,y,t]. If the function
    ! should be evaluated for a set of variables then ValBlock=[x,y] and
    ! ValScalar=[t] works fine.
!</description>

!<input>
    ! Function parser
    type (t_fparser),  intent(IN) :: rparser

    ! Function identifier
    integer, intent(IN) :: iComp

    ! Orientation of the stored values
    ! iDim =1 : ValBlock is organized as (x1:xN),(y1:yN),...
    ! iDim =2 : ValBlock is organized as (x1,y1),(x2,y2),...,(xN,yN)
    integer, intent(IN) :: iDim

    ! Variable values (must have the same dimension as Res)
    real(DP), dimension(:,:), intent(IN) :: ValBlock

    ! Variable values. This is a vector of scalar variables
    ! which is the same for all components of Res, e.g. the time variable.
    real(DP), dimension(:), intent(IN), optional :: ValScalar
!</input>

!<output>
    ! Evaluated function
    real(DP), dimension(:), intent(OUT) :: Res
!</output>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), allocatable :: ValTemp
    integer :: iVal,jVal,nVal,iMemory,iBlockSize,iBlock,sizeValBlock,sizeValScalar
    integer :: EvalErrType

    ! Get total number of variable sets
    nVal=size(ValBlock,iDim)

    ! Check if the compiled function is vectorizable
    if (rparser%Comp(iComp)%isVectorizable) then
      ! ...ok, vectorization of the bytecode is admissible.

      ! Guess required memory
      iMemory=(rparser%Comp(iComp)%StackSize+1)*nVal
      iMemory=min(iMemory,FPAR_MAXSTACKSIZE)
      
      ! Check if memory of global stack is sufficient for complete set of variables
      if (size(p_Stack) < iMemory) then
        if (rparser%Comp(iComp)%StackSize+1 < FPAR_MAXSTACKSIZE) then
          call storage_realloc('fparser_evalFunctionBlock',&
              iMemory,h_Stack,ST_NEWBLOCK_NOINIT,.false.)
          call storage_getbase_double(h_Stack,p_Stack)
        else
          call output_line('Stack size exceeds memory!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'fparser_evalFunctionBlock')
          call sys_halt()
        end if
      end if
      
      ! Compute size of blocks that can be vectorized
      iBlockSize = floor(real(iMemory)/real(rparser%Comp(iComp)%StackSize+1))
      
      ! What is the organization of ValBlock(:,:)
      if (iDim == 1) then
        do iVal=1,nVal,iBlockSize
          
          ! Initialization
          jVal     = min(iVal+iBlockSize,nVal)
          iBlock   = jVal-iVal+1
          
          ! Invoke working routine
          call evalFunctionBlock(iBlock,p_Stack,rparser%Comp(iComp),&
              ValBlock(iVal:jVal,:),iDim,EvalErrType,Res(iVal:jVal),ValScalar)
        end do
      else
        do iVal=1,nVal,iBlockSize
          
          ! Initialization
          jVal     = min(iVal+iBlockSize,nVal)
          iBlock   = jVal-iVal+1
          
          ! Invoke working routine
          call evalFunctionBlock(iBlock,p_Stack,rparser%Comp(iComp),&
              ValBlock(:,iVal:jVal),iDim,EvalErrType,Res(iVal:jVal),ValScalar)
        end do
      end if
      
    else   ! The compiled function cannot be vectorized

      ! Check if memory of global stack is sufficient for one set of variables
      if (size(p_Stack) < rparser%Comp(iComp)%StackSize+1) then
        if (rparser%Comp(iComp)%StackSize+1 < FPAR_MAXSTACKSIZE) then
          call storage_realloc('fparser_evalFunctionBlock',&
              rparser%Comp(iComp)%StackSize+1,h_Stack,ST_NEWBLOCK_NOINIT,.false.)
          call storage_getbase_double(h_Stack,p_Stack)
        else
          call output_line('Stack size exceeds memory!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'fparser_evalFunctionBlock')
          call sys_halt()
        end if
      end if

      ! The compiled bytecode cannot be vectorized. Hence, evaluate the function
      ! separately for each set of variables. Here, the organization of the array
      ! VAL(:,:) is important. Moreover, if the optional parameter ValScalar is
      ! given, then we have to combine those variables from ValBlock and ValScalar.

      if (present(ValScalar)) then

        ! Allocate auxiliary array
        sizeValBlock = size(ValBlock,3-iDim)
        sizeValScalar= size(ValScalar)
        allocate(ValTemp(sizeValBlock+sizeValScalar))

        if (iDim == 1) then
          do iVal=1,nVal
            
            ValTemp(1:sizeValBlock) = ValBlock(iVal,:)
            ValTemp(sizeValBlock+1:)= ValScalar

            ! Invoke working routine
            call evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValTemp,&
                EvalErrType,Res(iVal))
          end do
        else
          do iVal=1,nVal
            
            ValTemp(1:sizeValBlock) = ValBlock(:,iVal)
            ValTemp(sizeValBlock+1:)= ValScalar
            
            ! Invoke working routine
            call evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValTemp,&
                EvalErrType,Res(iVal))
          end do
        end if

        ! Deallocate auxiliary array
        deallocate(ValTemp)

      else
        
        if (iDim == 1) then
          do iVal=1,nVal
            
            ! Invoke working routine
            call evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValBlock(iVal,:),&
                EvalErrType,Res(iVal))
          end do
        else
          do iVal=1,nVal
            
            ! Invoke working routine
            call evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValBlock(:,iVal),&
                EvalErrType,Res(iVal))
          end do
        end if

      end if
    end if

    ! Check if evaluation was successful
    if (EvalErrType .ne. 0) then
      call output_line('An error occured during function evaluation!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'fparser_evalFunctionBlock')
      call sys_halt()
    end if
  end subroutine fparser_evalFunctionBlock

  ! *****************************************************************************

!<subroutine>

  subroutine evalFunctionScalar (Stack, Comp, Val, EvalErrType, res)

!<description>
    ! Evaluate bytecode for the values passed in array Val(:). Note, this subroutine
    ! requires some working memory Stack(*) which is not checked. So be warned,
    ! not to call this routine directly ;-)
!</description>

!<input>
    ! Component of function parser
    type(t_fparserComponent), intent(IN) :: Comp

    ! Variable values
    real(DP), dimension(:), intent(IN) :: Val
!</input>

!<inputoutput>
    ! Stack memory
    real(DP), dimension(*), intent(INOUT) :: Stack
!</inputoutput>

!<output>
    ! Error code for function evaluation
    integer, intent(OUT) :: EvalErrType

    ! Evaluated function
    real(DP), intent(OUT) :: res
!</output>
!</subroutine>

    ! local parameters
    real(DP), parameter :: Zero = 0._DP
    
    ! local variables
    integer  :: InstPtr,StackPtr,DataPtr
    integer  :: jumpAddr,immedAddr
    real(DP) :: daux
    
    ! Initialization
    DataPtr  = 1
    StackPtr = 0
    InstPtr  = 0
    
    ! Repeat until complete bytecode has been processed
    do while(InstPtr < Comp%ByteCodeSize)
      InstPtr=InstPtr+1

      ! What kind of bytecode are we?
      select case (Comp%ByteCode(InstPtr))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      case (cAbs)
        Stack(StackPtr)=abs(Stack(StackPtr))
        
      case (cAcos)
        if ((Stack(StackPtr) < -1._DP) .or. &
            (Stack(StackPtr) >  1._DP)) then
          EvalErrType = 4
          res         = zero
          return
        endif
        Stack(StackPtr)=acos(Stack(StackPtr))
        
      case (cAsin)
        if ((Stack(StackPtr) < -1._DP) .or. &
            (Stack(StackPtr) >  1._DP)) then
          EvalErrType = 4
          res         = zero
          return
        endif
        Stack(StackPtr)=asin(Stack(StackPtr))
        
      case (cAtan)
        Stack(StackPtr)=atan(Stack(StackPtr))

      case (cAtan2)
        Stack(StackPtr-1)=atan2(Stack(StackPtr -1),Stack(StackPtr))
        StackPtr=StackPtr-1
        
      case (cAcosh)
        daux=Stack(StackPtr)+sqrt(Stack(StackPtr)**2-1)
        if (daux <= 0._DP) then
          EvalErrType = 5
          res         = zero
          return
        end if
        Stack(StackPtr)=log(daux)
        
      case (cAnint)
        Stack(StackPtr)=anint(Stack(StackPtr))

      case (cAint)
        Stack(StackPtr)=aint(Stack(StackPtr))

      case (cAsinh)
        daux=Stack(StackPtr)+sqrt(Stack(StackPtr)**2+1)
        if (daux <= 0._DP) then
          EvalErrType = 5
          res         = zero
          return
        end if
        Stack(StackPtr)=log(daux)
        
      case (cAtanh)
        if (Stack(StackPtr) == -1._DP) then
          EvalErrType = 6
          res         = zero
          return
        end if
        daux=(1+Stack(StackPtr))/(1-Stack(StackPtr))
        if (daux <= 0._DP) then
          EvalErrType = 3
          res         = zero
          return
        end if
        Stack(StackPtr)=log(daux)/2._DP
        
      case (cCeil)
        Stack(StackPtr)=ceiling(Stack(StackPtr))
        
      case (cCos)
        Stack(StackPtr)=cos(Stack(StackPtr)) 
        
      case (cCosh)
        Stack(StackPtr)=cosh(Stack(StackPtr))

      case (cCot)
        daux=tan(Stack(StackPtr))
        if (daux == 0._DP) then 
          EvalErrType = 1
          res         = zero
          return
        end if
        Stack(StackPtr)=1._DP/daux

      case (cCsc)
        daux=sin(Stack(StackPtr))
        if (daux == 0._DP) then
          EvalErrType = 1
          res         = zero
          return
        endif
        Stack(StackPtr)=1._DP/daux
        
      case (cExp)
        Stack(StackPtr)=exp(Stack(StackPtr))

      case (cFloor)
        Stack(StackPtr)=floor(Stack(StackPtr))

      case (cIf)
        InstPtr=InstPtr+1; jumpAddr  = Comp%ByteCode(InstPtr)
        InstPtr=InstPtr+1; immedAddr = Comp%ByteCode(InstPtr)
        if (.not.DbleToLogc(Stack(StackPtr))) then
          InstPtr = jumpAddr
          DataPtr = immedAddr
        end if
        StackPtr=StackPtr-1
        
      case (cLog)
        if (Stack(StackPtr) <= 0._DP) then
          EvalErrType = 3
          res         = zero
          return
        endif
        Stack(StackPtr)=log(Stack(StackPtr)) 
        
      case (cLog10)
        if (Stack(StackPtr) <= 0._DP) then
          EvalErrType = 3
          res         = zero
          return
        endif
        Stack(StackPtr)=log10(Stack(StackPtr))
        
      case (cMax)
        Stack(StackPtr-1)=max(Stack(StackPtr-1),Stack(StackPtr))
        StackPtr=StackPtr-1
        
      case (cMin)
        Stack(StackPtr-1)=min(Stack(StackPtr-1),Stack(StackPtr))
        StackPtr=StackPtr-1
        
      case (cSec)
        daux=cos(Stack(StackPtr))
        if (daux == 0._DP) then
          EvalErrType = 1
          res         = zero
          return
        endif
        Stack(StackPtr)=1._DP/daux
        
      case (cSin)
        Stack(StackPtr)=sin(Stack(StackPtr))
        
      case(cSinh)
        Stack(StackPtr)=sinh(Stack(StackPtr))
      
      case(cSqrt)
        if (Stack(StackPtr) < 0._DP) then
          EvalErrType = 3
          res         = zero
          return
        endif
        Stack(StackPtr)=sqrt(Stack(StackPtr))

      case (cTan)
        Stack(StackPtr)=tan(Stack(StackPtr))
        
      case (cTanh)
        Stack(StackPtr)=tanh(Stack(StackPtr))
        
        !------------------------------------------------------------
        ! Misc
        !------------------------------------------------------------
      case (cImmed)
        StackPtr        = StackPtr+1
        Stack(StackPtr) = Comp%Immed(DataPtr)
        DataPtr         = DataPtr+1
        
      case (cJump)
        DataPtr=Comp%ByteCode(InstPtr+2)
        InstPtr=Comp%ByteCode(InstPtr+1)

        !------------------------------------------------------------
        ! Operators
        !------------------------------------------------------------
      case (cNeg)
        Stack(StackPtr)=-Stack(StackPtr)
        
      case (cAdd)
        Stack(StackPtr-1)=Stack(StackPtr-1)+Stack(StackPtr)
        StackPtr=StackPtr-1
        
      case (cSub)
        Stack(StackPtr-1)=Stack(StackPtr-1)-Stack(StackPtr)
        StackPtr=StackPtr-1

      case (cMul)
        Stack(StackPtr-1)=Stack(StackPtr-1)*Stack(StackPtr)
        StackPtr=StackPtr-1
        
      case (cDiv)
        if (Stack(StackPtr) == 0._DP) then
          EvalErrType = 1
          res         = zero
          return
        endif
        Stack(StackPtr-1)=Stack(StackPtr-1)/Stack(StackPtr)
        StackPtr=StackPtr-1
        
      case (cMod)
        if (Stack(StackPtr) == 0._DP) then
          EvalErrType = 1
          res         = zero
          return
        endif
        Stack(StackPtr-1)=mod(Stack(StackPtr-1),Stack(StackPtr))
        StackPtr=StackPtr-1

      case (cPow)
        Stack(StackPtr-1)=Stack(StackPtr-1)**Stack(StackPtr)
        StackPtr=StackPtr-1
        
      case (cEqual)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) == Stack(StackPtr))
        StackPtr=StackPtr-1

      case (cNEqual)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) /= Stack(StackPtr))
        StackPtr=StackPtr-1

      case (cLess)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) < Stack(StackPtr))
        StackPtr=StackPtr-1

      case (cLessOrEq)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) <= Stack(StackPtr))
        StackPtr=StackPtr-1
        
      case (cGreater)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) > Stack(StackPtr))
        StackPtr=StackPtr-1
        
      case (cGreaterOrEq)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) >= Stack(StackPtr))
        StackPtr=StackPtr-1
        
      case (cAnd)
        Stack(StackPtr-1)=LogcToDble(DbleToLogc( Stack(StackPtr-1)) .and. &
            DbleToLogc(Stack(StackPtr)) )
        StackPtr=StackPtr-1

      case (cOr)
        Stack(StackPtr-1)=LogcToDble(DbleToLogc( Stack(StackPtr-1)) .or. &
            DbleToLogc(Stack(StackPtr)) )
        StackPtr=StackPtr-1

      case (cNot)
        Stack(StackPtr)=LogcToDble( .not. DbleToLogc(Stack(StackPtr)) )
        
        !------------------------------------------------------------
        ! Degrees-radians conversion
        !------------------------------------------------------------
      case (cDeg)
        Stack(StackPtr)=RadToDeg(Stack(StackPtr))
        
      case (cRad)
        Stack(StackPtr)=DegToRad(Stack(StackPtr))
        
      case DEFAULT
        StackPtr=StackPtr+1
        Stack(StackPtr)=Val(Comp%ByteCode(InstPtr)-VarBegin+1)
      end select
    end do
    
    EvalErrType = 0
    Res = Stack(StackPtr)
  end subroutine evalFunctionScalar

  ! *****************************************************************************

!<subroutine>

  subroutine evalFunctionBlock (BlockSize, Stack, Comp, ValBlock, iDim, EvalErrType, Res, ValScalar)

!<description>
    ! Evaluate bytecode for an array of values passed in ValBlock(:,:).
    ! Note, this subroutine requires some working memory Stack(iBlock,*) which is
    ! not checked. So be warned, not to call this function directly ;-)
!</description>

!<input>
    ! Component of function parser
    type(t_fparserComponent), intent(IN) :: Comp
    
    ! Variable values
    real(DP), dimension(:,:), intent(IN) :: ValBlock

    ! Size of the vector block
    integer, intent(IN) :: BlockSize

    ! Orientation of the stored values
    ! iDim =1 : ValBlock is organized as (x1:xN),(y1:yN),...
    ! iDim =2 : ValBlock is organized as (x1,y1),(x2,y2),...,(xN,yN)
    integer, intent(IN) :: iDim

    ! Vector of scalar variable values
    real(DP), dimension(:), intent(IN), optional :: ValScalar
!</input>

!<inputoutput>
    ! Stack memory
    real(DP), dimension(BlockSize,*), intent(INOUT) :: Stack
!</inputoutput>

!<output>
    ! Error code for function evaluation
    integer, intent(OUT) :: EvalErrType

    ! Evaluated function
    real(DP), dimension(:) :: Res
!</output>
!</subroutine>

    ! local parameters
    real(DP), parameter :: Zero = 0._DP
    
    ! local variables
    integer  :: InstPtr,DataPtr,StackPtr,iBlock,istartValScalar,iVariable
    real(DP) :: daux
    
    ! Initialization
    DataPtr  = 1
    StackPtr = 0
    InstPtr  = 0

    ! This is tricky. istartValScalar indicates the number of the first
    ! variable which is passed as scalar. Hence, if the optional parameter
    ! ValScalar is missing, then istartValScalar pointers to SIZE(ValBlock)+1.
    ! Obviously, no variable beyond this value is addressed.
    istartValScalar=size(ValBlock,3-iDim)+1
    
    ! Repeat until complete bytecode has been processed
    do while(InstPtr < Comp%ByteCodeSize)
      InstPtr=InstPtr+1
      
      select case (Comp%ByteCode(InstPtr))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      case (cAbs)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=abs(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        

      case (cAcos)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
          if ((Stack(iBlock,StackPtr) < -1._DP) .or.&
              (Stack(iBlock,StackPtr) >  1._DP)) then
            EvalErrType=4
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr)=acos(Stack(iBlock,StackPtr))
          end if
        end do
!$omp end parallel do


      case (cAsin)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
          if ((Stack(iBlock,StackPtr) < -1._DP) .or.&
              (Stack(iBlock,StackPtr) >  1._DP)) then
            EvalErrType=4
            Res(iBlock)=zero
          else
            Stack(iBlock,StackPtr)=asin(Stack(iBlock,StackPtr))
          end if
        end do
!$omp end parallel do

        
      case (cAtan)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=atan(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do


      case (cAtan2)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=atan2(Stack(iBlock,StackPtr -1),Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1

        
      case (cAcosh)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
          daux=Stack(iBlock,StackPtr)+sqrt(Stack(iBlock,StackPtr)**2-1)
          if (daux <= 0._DP) then
            EvalErrType=5
            Res(iBlock)=zero
          else
            Stack(iBlock,StackPtr)=log(daux)
          end if
        end do
!$omp end parallel do
        

      case (cAnint)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=anint(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do


      case (cAint)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=aint(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do


      case (cAsinh)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
          daux=Stack(iBlock,StackPtr)+sqrt(Stack(iBlock,StackPtr)**2+1)
          if (daux <= 0._DP) then
            EvalErrType=5
            Res(iBlock)=zero
          else
            Stack(iBlock,StackPtr)=log(daux)
          end if
        end do
!$omp end parallel do


      case (cAtanh)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
          if (Stack(iBlock,StackPtr) == -1._DP) then
            EvalErrType=6
            Res(iBlock)=Zero
          end if
          daux=(1+Stack(iBlock,StackPtr))/(1-Stack(iBlock,StackPtr))
          if (daux <= 0._DP) then
            EvalErrType=3
            Res(iBlock)=zero
          else
            Stack(iBlock,StackPtr) = log(daux)/2._DP
          end if
        end do
!$omp end parallel do

        
      case (cCeil)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=ceiling(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do


      case (cCos)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=cos(Stack(iBlock,StackPtr)) 
        end do
!$omp end parallel do
        

      case (cCosh)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=cosh(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do


      case (cCot)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
          daux=tan(Stack(iBlock,StackPtr))
          if (daux == 0) then 
            EvalErrType=1
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr)=1._DP/daux
          end if
        end do
!$omp end parallel do


      case (cCsc)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
          daux=sin(Stack(iBlock,StackPtr))
          if (daux==0._DP) then
            EvalErrType=1
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr)=1._DP/daux
          end if
        end do
!$omp end parallel do


      case (cExp)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=exp(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do


      case (cFloor)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=floor(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do


      case (cIf)
        ! IF-THEN-ELSE cannot be vectorized which should be noted during
        ! bytecode compilation. If we reach this point, then something
        ! went wrong before.
        call output_line('IF-THEN-ELSE cannot be vectorized!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'evalFunctionBlock')
        call sys_halt()
        

      case (cLog)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
          if (Stack(iBlock,StackPtr) <= 0._DP) then
            EvalErrType=3
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr)=log(Stack(iBlock,StackPtr)) 
          end if
        end do
!$omp end parallel do


      case (cLog10)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
          if (Stack(iBlock,StackPtr) <= 0._DP) then
            EvalErrType=3
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr)=log10(Stack(iBlock,StackPtr))
          end if
        end do
!$omp end parallel do


      case (cMax)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=max(Stack(iBlock,StackPtr-1),Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1
        

      case (cMin)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=min(Stack(iBlock,StackPtr-1),Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1
        

      case (cSec)
!$omp parallel do default(shared) private(iBlock,daux)
        do iBlock=1,BlockSize
          daux=cos(Stack(iBlock,StackPtr))
          if (daux == 0._DP) then
            EvalErrType=1
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr)=1._DP/daux
          end if
        end do
!$omp end parallel do

        
      case (cSin)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=sin(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do

        
      case(cSinh)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=sinh(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
      

      case(cSqrt)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
          if (Stack(iBlock,StackPtr) < 0._DP) then
            EvalErrType=3
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr)=sqrt(Stack(iBlock,StackPtr))
          end if
        end do
!$omp end parallel do


      case (cTan)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=tan(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do

        
      case (cTanh)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=tanh(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        

        !------------------------------------------------------------
        ! Misc
        !------------------------------------------------------------
      case (cImmed)
        StackPtr=StackPtr+1
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=Comp%Immed(DataPtr)
        end do
!$omp end parallel do
        DataPtr=DataPtr+1


      case (cJump)
        DataPtr=Comp%ByteCode(InstPtr+2)
        InstPtr=Comp%ByteCode(InstPtr+1)


        !------------------------------------------------------------
        ! Operators
        !------------------------------------------------------------
      case (cNeg)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=-Stack(iBlock,StackPtr)
        end do
!$omp end parallel do


      case (cAdd)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)+Stack(iBlock,StackPtr)
        end do
!$omp end parallel do
        StackPtr=StackPtr-1
        

      case (cSub)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)-Stack(iBlock,StackPtr)
        end do
!$omp end parallel do
        StackPtr=StackPtr-1


      case (cMul)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)*Stack(iBlock,StackPtr)
        end do
!$omp end parallel do
        StackPtr=StackPtr-1
        

      case (cDiv)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
          if (Stack(iBlock,StackPtr) == 0._DP) then
            EvalErrType=1
            Res(iBlock)=zero
          else
            Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)/Stack(iBlock,StackPtr)
          end if
        end do
!$omp end parallel do
        StackPtr=StackPtr-1
        

      case (cMod)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
          if (Stack(iBlock,StackPtr) == 0._DP) then
            EvalErrType=1
            Res(iBlock)=Zero
          else
            Stack(iBlock,StackPtr-1)=mod(Stack(iBlock,StackPtr-1),Stack(iBlock,StackPtr))
          end if
        end do
!$omp end parallel do
        StackPtr=StackPtr-1


      case (cPow)
!$omp parallel do default(shared) private(iBlock)
        do  iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)**Stack(iBlock,StackPtr)
        end do
!$omp end parallel do
        StackPtr=StackPtr-1

        
      case (cEqual)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) == Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1


      case (cNEqual)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) /= Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1


      case (cLess)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) < Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1


      case (cLessOrEq)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) <= Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1

        
      case (cGreater)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) > Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1

        
      case (cGreaterOrEq)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) >= Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        StackPtr=StackPtr-1

        
      case (cAnd)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(DbleToLogc(Stack(iBlock,StackPtr-1)) .and. &
                DbleToLogc(Stack(iBlock,StackPtr)) )
        end do
!$omp end parallel do
        StackPtr=StackPtr-1


      case (cOr)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(DbleToLogc(Stack(iBlock,StackPtr-1)) .or. &
                DbleToLogc(Stack(iBlock,StackPtr)) )
        end do
!$omp end parallel do
        StackPtr=StackPtr-1


      case (cNot)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=LogcToDble( .not. DbleToLogc(Stack(iBlock,StackPtr)) )
        end do
!$omp end parallel do

        
        !------------------------------------------------------------
        ! Degrees-radians conversion
        !------------------------------------------------------------
      case (cDeg)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=RadToDeg(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do

        
      case (cRad)
!$omp parallel do default(shared) private(iBlock)
        do iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=DegToRad(Stack(iBlock,StackPtr))
        end do
!$omp end parallel do
        

      case DEFAULT
        StackPtr  = StackPtr+1
        iVariable = Comp%ByteCode(InstPtr)-VarBegin+1

        ! Do we have to process one of the scalar variables of one of the block variables
        if (iVariable >= istartValScalar) then
!$omp parallel do default(shared) private(iBlock)
          do iBlock=1,BlockSize
             Stack(iBlock,StackPtr)=ValScalar(iVariable-istartValScalar+1)
          end do
!$omp end parallel do
        else
          if (iDim == 1) then
!$omp parallel do default(shared) private(iBlock)
            do iBlock=1,BlockSize
               Stack(iBlock,StackPtr)=ValBlock(iBlock,iVariable)
            end do
!$omp end parallel do
          else
!$omp parallel do default(shared) private(iBlock)
            do iBlock=1,BlockSize
               Stack(iBlock,StackPtr)=ValBlock(iVariable,iBlock)
            end do
!$omp end parallel do
          end if
        end if
      end select
    end do
    
    EvalErrType = 0
!$omp parallel do default(shared) private(iBlock)
    do iBlock=1,BlockSize
       Res(iBlock) = Stack(iBlock,StackPtr)
    end do
!$omp end parallel do
  end subroutine evalFunctionBlock
  
  ! *****************************************************************************

!<function>

  function fparser_ErrorMsg (EvalErrType) result (msg)

!<description>
    ! Return error message of function parser
!</description>

    ! local constants
    character (LEN=*), dimension(6), parameter :: m = (/ 'Division by zero                  ', &
                                                         'Argument of SQRT negative         ', &
                                                         'Argument of LOG negative          ', &
                                                         'Argument of ASIN or ACOS illegal  ', &
                                                         'Argument of ASINH or ACOSH illegal', &
                                                         'Argument of ATANH illegal         ' /)

!<input>
    ! Error identifier
    integer, intent(IN) :: EvalErrType
!</input>   
    
!<result>
    ! Error messages
    character (LEN=len(m)) :: msg
!</result>
!</function>

    if (EvalErrType < 1 .or. EvalErrType > size(m)) then
      msg = ''
    else
      msg = m(EvalErrType)
    endif
  end function fparser_ErrorMsg

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_PrintByteCode(rparser, iComp)

!<description>
    ! Print the compiled bytecode stack
!</description>

!<input>
    ! Function parser
    type(t_fparser), intent(IN) :: rparser

    ! Function identifier
    integer, intent(IN) :: iComp
!</input>
!</subroutine>

    ! local variables
    type(t_fparserComponent), pointer :: Comp
    character(LEN=5) :: n
    integer :: InstPtr, DataPtr, StackPtr,nparams
    integer(is) :: opCode
    
    nparams  = 1
    DataPtr  = 1
    StackPtr = 0
    InstPtr  = 0
    Comp     => rparser%Comp(iComp)

    do while(InstPtr < Comp%ByteCodeSize)
      InstPtr=InstPtr+1
      
      write(*,FMT='(I8.8,1X,":",1X)',ADVANCE="NO") InstPtr
      
      opCode=Comp%ByteCode(InstPtr)
      select case(opCode)
        
      case(cIf)
        write(*,FMT='(A,1X,T10,I8.8)') "jz",Comp%ByteCode(InstPtr+1)+1
        InstPtr=InstPtr+2

      case(cJump)
        write(*,FMT='(A,1X,T10,I8.8)') "jump",Comp%ByteCode(InstPtr+1)+1
        InstPtr=InstPtr+2

      case(cImmed)
        write(*,FMT='(A,1X,T10,G8.2)') "push",Comp%Immed(DataPtr)
        DataPtr=DataPtr+1

      case DEFAULT
        if (opCode < VarBegin) then
          select case(opCode)
          case(cNEG); n="neg"
          case(cADD); n="add"
          case(cSUB); n="sub"
          case(cMUL); n="mul"
          case(cDIV); n="div"
          case(cMOD); n="mod"
          case(cPOW); n="pow"
          case(cEqual); n="eq"
          case(cNEqual); n="ne"
          case(cLess); n="lt"
          case(cLessOrEq); n="le"
          case(cGreater); n="gt"
          case(cGreaterOrEq); n="ge"
          case(cAND); n="and"
          case (cOR); n="or"
          case(cNOT); n="not"
          case(cDEG); n="deg"
          case(cRAD); n="rad"

          case DEFAULT
            n       = Funcs(opCode)
            nparams = MathFunctionParameters(opCode)
          end select
          write(*,FMT='(A,T10,A,"  (",I1,") ")') trim(n),"Par",nparams
          
        else
          write(*,FMT='(A,T10,A,1X,I4.4)') "push","Var",opCode-VarBegin+1
        end if
        
      end select
    end do
  end subroutine fparser_PrintByteCode

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<subroutine>
  
  recursive subroutine CheckSyntax (Func,Var)

!<description>
    ! Check syntax of function string, returns 0 if syntax is ok
!</description>
    
!<input>
    ! Function string without spaces
    character (LEN=*), intent(in) :: Func

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Var
!</input>
!</subroutine>
    
    ! local avriables
    type(t_stack) :: functionParenthDepth
    integer(is) :: n,opSize
    character (LEN=1) :: c
    real(DP) :: r
    logical :: err
    integer :: FuncIndex,FuncIndex2
    integer :: ParCnt,ib,in,FuncLength,idummy
    
    ! Initialization
    FuncIndex = 1
    ParCnt = 0
    FuncLength = len_trim(Func)
    call stack_create(functionParenthDepth,max(5,int(FuncLength/4._DP)),ST_INT)

    do
      if (FuncIndex > FuncLength) then
        call output_line('Invalid function string!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
        call sys_halt()
      end if
      c = Func(FuncIndex:FuncIndex)

      ! Check for valid operand (must appear)

      ! Check for leading - or !
      if (c == '-' .or. c == '!') then                      
        FuncIndex = FuncIndex+1
        if (FuncIndex > FuncLength) then
          call output_line('Premature end of string!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        c = Func(FuncIndex:FuncIndex)
      end if
            
      ! Check for math function
      n = MathFunctionIndex (Func(FuncIndex:))
      if (n > 0) then
        ! Math function found
        FuncIndex = FuncIndex+len_trim(Funcs(n))
        if (FuncIndex > FuncLength) then
          call output_line('Premature end of string!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        c = Func(FuncIndex:FuncIndex)
        if (c /= '(') then
          call output_line('Expecting ( after function '//Func(FuncIndex:)//'!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        FuncIndex2=FuncIndex+1
        if (FuncIndex2 > FuncLength) then
          call output_line('Premature end of string!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        if (Func(FuncIndex2:FuncIndex2) == ')') then
          FuncIndex=FuncIndex2+1
          if (FuncIndex > FuncLength) then
            call output_line('Premature end of string!',&
                             OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
            call sys_halt()
          end if
          c = Func(FuncIndex:FuncIndex)
          
          ! Ugly, but other methods would just be uglier ...
          goto 999
        end if
        
        ! Push counter for parenthesss to stack
        call stack_pushback(functionParenthDepth,ParCnt+1)
      end if
      
      ! Check for opening parenthesis
      if (c == '(') then
        ParCnt = ParCnt+1
        FuncIndex = FuncIndex+1
        if (FuncIndex > FuncLength) then
          call output_line('Premature end of string!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        if (Func(FuncIndex:FuncIndex) == ')') then
          call output_line('Empty parantheses!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        cycle
      end if

      ! Check for number
      if (scan(c,'0123456789.') > 0) then
        r = RealNum (Func(FuncIndex:),ib,in,err)
        if (err) then
          call output_line('Invalid number format!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        FuncIndex = FuncIndex+in-1
        if (FuncIndex > FuncLength) exit
        c = Func(FuncIndex:FuncIndex)
      elseif (c == '_') then
        ! Check for constant
        n = ConstantIndex (Func(FuncIndex:))
        if (n == 0) then
          call output_line('Invalid constant!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        FuncIndex = FuncIndex+len_trim(Consts(n))+1
        if (FuncIndex > FuncLength) exit
        c=Func(FuncIndex:FuncIndex)
      elseif (c == '@') then
        ! Check for expression
        n = ExpressionIndex (Func(FuncIndex:))
        if (n == 0) then
          call output_line('Invalid expression!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        call CheckSyntax(ExpressionStrings(n), Var)
        FuncIndex = FuncIndex+len_trim(Expressions(n))+1
        if (FuncIndex > FuncLength) exit
        c=Func(FuncIndex:FuncIndex)
      else
        ! Check for variable
        n = VariableIndex (Func(FuncIndex:),Var,ib,in)
        if (n == 0) then
          call output_line('Invalid element!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        FuncIndex = FuncIndex+in-1
        if (FuncIndex > FuncLength) exit
        c = Func(FuncIndex:FuncIndex)
      end if

      ! Check for closing parenthesis
      do while (c == ')')
        if (.not.stack_isempty(functionParenthDepth)) then
          if(stack_backInt(functionParenthDepth) == ParCnt) &
              idummy=stack_popbackInt(functionParenthDepth)
        end if
        ParCnt = ParCnt-1
        if (ParCnt < 0) then
          call output_line('Mismatched parenthesis!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        if (Func(FuncIndex-1:FuncIndex-1) == '(') then
          call output_line('Empty parentheses!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        FuncIndex = FuncIndex+1
        if (FuncIndex > FuncLength) exit
        c = Func(FuncIndex:FuncIndex)
      end do

      ! Now, we have a legal operand: A legal operator or end of string must follow
999   if (FuncIndex > FuncLength) exit
      
      ! Check operators
      opSize = 0
      if (.not.stack_isempty(functionParenthDepth)) then
        if (c == ',' .and. &
            stack_backInt(functionParenthDepth) == parCnt) then
          opSize = 1
        else
          opSize = isOperator(Func(FuncIndex:))
        end if
      else
        opSize = isOperator(Func(FuncIndex:))
      end if
      if (opSize == 0) then
        call output_line('Operator expected!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
        call sys_halt()
      end if
      
      ! Now, we have an operand and an operator: the next loop will check for another 
      ! operand (must appear)
      FuncIndex = FuncIndex+opSize
    end do
    if (ParCnt > 0) then
      call output_line('Missing )!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'CheckSyntax')
      call sys_halt()
    end if
    call stack_release(functionParenthDepth)
  end subroutine CheckSyntax

  ! *****************************************************************************

!<function>

  function isOperator (str) result (n)

!<description>
    ! Return 0 if given string is not an operator, else the size of the
    ! operator
!</description>
    
!<input>
    ! Operator string
    character(LEN=*), intent(IN) :: str    
!</input>

!<result>
    ! Length of operator, 0 if string is no operator
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j,m

    n = 0
    do j=cAdd,cOr
      m = len_trim(Ops(j))
      if (str(1:m) == trim(Ops(j))) then
        n = m
        exit
      end if
    end do
    
  end function isOperator

  ! *****************************************************************************

!<function>

  function MathFunctionIndex (str) result (n)

!<description>
    ! Return index of math function beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    character (LEN=*), intent(in) :: str
!</input>

!<result>
    ! Index of math function
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j
    integer k
    character (LEN=len(Funcs)) :: fun
    
    ! Check all math functions
    n = 0
    do j=cIf,cSec
      k = min(len_trim(Funcs(j)), len(str))

      ! Compare lower case letters
      call sys_tolower(str(1:k), fun)
      if (fun == Funcs(j)) then
        ! Found a matching function
        n = j
        return
      end if
    end do
  end function MathFunctionIndex

  ! *****************************************************************************
  
!<function>

  function MathFunctionParameters (FuncIndex) result (nparameters)

!<description>
    ! Return number of required parameters
!</description>

!<input>
    ! Index of function
    integer(is) :: FuncIndex
!</input>

!<result>
    ! Number if required parameters
    integer :: nparameters
!</result>
!</function>

    select case(FuncIndex)
    case(cIf)
      nparameters = 3
    case(cMin:cAtan2)
      nparameters = 2
    case(cAbs:cSec)
      nparameters = 1
    case DEFAULT
      nparameters=0
      call output_line('Not a function',&
                       OU_CLASS_WARNING,OU_MODE_STD,'MathFunctionParameters')
    end select
  end function MathFunctionParameters

  ! *****************************************************************************

!<function>

  function ConstantIndex (str) result (n)

!<description>
    ! Return index of predefined constants beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    character (LEN=*), intent(in) :: str
!</input>

!<result>
    ! Index of math function
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j
    integer k
    character (LEN=len(Consts)) :: con
    
    ! Check all predefined constants
    n = 0
    do j=1,nConsts
      k = min(len_trim(Consts(j)), len(str(2:)))

      ! Compare lower case letters
      call sys_tolower(str(2:k+1),con)
      if (con == Consts(j)) then
        ! Found a matching constant
        n = j
        return
      end if
    end do
  end function ConstantIndex

  ! *****************************************************************************

!<function>

  function ExpressionIndex (str) result (n)

!<description>
    ! Return index of predefined expression beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    character (LEN=*), intent(in) :: str
!</input>

!<result>
    ! Index of math function
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j
    integer k
    character (LEN=len(Expressions)) :: expression

    ! Check all predefined expressions
    n = 0
    do j=1,nExpressions
      k = min(len_trim(Expressions(j)), len(str(2:)))

      ! Compare lower case letters
      call sys_tolower(str(2:k+1),expression)

      if (expression == Expressions(j)) then
        ! Found a matching expression
        n = j
        return
      end if
    end do
  end function ExpressionIndex

  ! *****************************************************************************

!<function>
  
  function VariableIndex (str, Var, ibegin, inext) result (n)

!<description>
    ! Return index of variable at begin of string str (returns 0 if no variable found)
!</description>

!<input>
    ! String
    character (LEN=*), intent(in) :: str

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Var
!</input>

!<output>
    ! OPTIONAL: Start position of variable name
    integer, optional, intent(out) :: ibegin

    ! OPTIONAL: Position of character after name
    integer, optional, intent(out) :: inext
!</output>

!<result>
    ! Index of variable
    integer(is) :: n
!</result>
!</function>

    ! local variables    
    integer :: j,ib,in,lstr
    
    n = 0
    lstr = len_trim(str)
    if (lstr > 0) then
      ! Search for first character in str
      do ib=1,lstr
        ! When lstr>0 at least 1 char in str
        if (str(ib:ib) /= ' ') exit
      end do
      
      ! Search for name terminators
      do in=ib,lstr
        if (scan(str(in:in),'+-*/%^),&|<>=! ') > 0) exit
      end do
      do j=1,size(Var)
        if (str(ib:in-1) == Var(j)) then
          ! Variable name found
          n = j
          exit
        end if
      end do
    end if

    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
  end function VariableIndex

  ! *****************************************************************************

!<subroutine>

  subroutine RemoveSpaces (str)

!<description>
    ! Remove Spaces from string, remember positions of characters in
    ! old string
!</description>
   
!<inputoutput>
    ! String from which spaces should be removed
    character (LEN=*), intent(inout) :: str
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: k,lstr
    
    lstr = len_trim(str)
    k = 1
    do while (str(k:lstr) /= ' ')                             
      if (str(k:k) == ' ') then
        str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
        k = k-1
      end if
      k = k+1
    end do
  end subroutine RemoveSpaces

  ! *****************************************************************************

!<subroutine>

  subroutine Replace (ca,cb,str)

!<description>
    ! Replace ALL appearances of character set ca in string str by character set cb
!</description>

!<input>
    ! Source characters
    character (LEN=*), intent(in) :: ca
    
    ! Destination characters
    character (LEN=len(ca)), intent(in) :: cb
!</input>

!<inputoutput>
    ! String
    character (LEN=*), intent(inout) :: str
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: j,lca
    
    lca = len(ca)
    do j=1,len_trim(str)-lca+1
      if (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    end do
  end subroutine Replace

  ! *****************************************************************************

!<subroutine>

  subroutine Compile (Comp, Func, Var)

!<description>
    ! Compile i-th function string Func into bytecode
!</description>

!<input>
    ! Function string
    character (LEN=*), intent(in) :: Func

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Var
!</input>

!<inputoutput>
    ! Function parser component
    type (t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    ! local variables
    integer(is), dimension(:), pointer :: ByteCode
    real(DP), dimension(:), pointer :: Immed
    integer :: ind,isize
    
    ! (Re-)initialize the bytecode structure (if required)
    if (associated(Comp%ByteCode)) deallocate (Comp%ByteCode)
    if (associated(Comp%Immed))    deallocate (Comp%Immed)
    Comp%ByteCodeSize = 0
    Comp%ImmedSize    = 0
    Comp%StackSize    = 0
    Comp%StackPtr     = 0
    Comp%isVectorizable = .true.

    ! Neither the stack for the bytecode nor the stack for the immediate expressions
    ! can exceed the size of the function string. Hence, allocate some initial memory
    isize=FunctionSize(Func)
    allocate(Comp%ByteCode(isize),Comp%Immed(isize))

    ! Compile function string into bytecode
    ind = CompileExpression(Comp,Func,1,Var)
    
    ! Adjust memory size of bytecode stack
    allocate(ByteCode(Comp%ByteCodeSize))
    ByteCode=Comp%ByteCode
    deallocate(Comp%ByteCode)
    allocate(Comp%ByteCode(Comp%ByteCodeSize))
    Comp%ByteCode=ByteCode
    deallocate(ByteCode)
    
    ! Adjust memory size of immediate stack
    allocate(Immed(Comp%ImmedSize))
    Immed=Comp%Immed
    deallocate(Comp%Immed)
    allocate(Comp%Immed(Comp%ImmedSize))
    Comp%Immed=Immed
    deallocate(Immed)
  end subroutine Compile

  ! *****************************************************************************

!<subroutine>
  
  subroutine incStackPtr (Comp)

!<description>
    ! Increase stack pointer
!</description>

!<inputoutput>
    type (t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    Comp%StackPtr=Comp%StackPtr+1
    if (Comp%StackPtr > Comp%StackSize) Comp%StackSize=Comp%StackSize+1
  end subroutine incStackPtr

  ! *****************************************************************************

!<subroutine>

  subroutine AddCompiledByte (Comp, byte)

!<description>
    ! Add compiled byte to bytecode
!</description>

!<input>
    ! Value of byte to be added
    integer(is), intent(in) :: byte
!</input>

!<inputoutput>
    type (t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP) :: daux

    Comp%ByteCodeSize = Comp%ByteCodeSize + 1    
    Comp%ByteCode(Comp%ByteCodeSize) = byte

    ! Try to optimize the compiled bytecode. Check the bytecode instruction and
    ! compute some values on-the-fly of this is possible
    select case(byte)
      !------------------------------------------------------------
      ! Functions
      !------------------------------------------------------------
    case (cAbs)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=abs(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cAcos)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) < -1._DP .or. &
            Comp%Immed(Comp%ImmedSize) >  1._DP) then
          call output_line('Invalid argument for ACOS!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=acos(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if
      
    case (cAsin)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) < -1._DP .or. &
            Comp%Immed(Comp%ImmedSize) >  1._DP) then
          call output_line('Invalid argument for ASIN!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=asin(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if
      
    case (cAtan)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=atan(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cAtan2)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=atan2(Comp%Immed(Comp%ImmedSize),Comp%Immed(Comp%ImmedSize-1))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cAcosh)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        daux=Comp%Immed(Comp%ImmedSize)+sqrt(Comp%Immed(Comp%ImmedSize)**2-1)
        if (daux <= 0) then
          call output_line('Invalid argument for ACOSH!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=log(daux)
        call RemoveCompiledByte(Comp)
      end if
      
    case (cAnint)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=anint(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cAint)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=aint(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cAsinh)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        daux=Comp%Immed(Comp%ImmedSize)+sqrt(Comp%Immed(Comp%ImmedSize)**2-1)
        if (daux <= 0) then
          call output_line('Invalid argument for ASINH!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=log(daux)
        call RemoveCompiledByte(Comp)
      end if

    case (cAtanh)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) == -1._DP) then
          call output_line('Invalid argument for ATANH!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        daux=(1+Comp%Immed(Comp%ImmedSize))/(1-Comp%Immed(Comp%ImmedSize))
        if (daux <= 0._DP) then
          call output_line('Invalid argument for ATANH!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=log(daux)/2._DP
        call RemoveCompiledByte(Comp)
      end if

    case (cCeil)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=ceiling(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cCos)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=cos(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cCosh)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=cosh(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cCot)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        daux=tan(Comp%Immed(Comp%ImmedSize))
        if (daux == 0._DP) then
          call output_line('Invalid argument for COT!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=1/daux
        call RemoveCompiledByte(Comp)
      end if
      
    case (cCsc)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        daux=sin(Comp%Immed(Comp%ImmedSize))
        if (daux == 0._DP) then
          call output_line('Invalid argument for CSC!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=1/daux
        call RemoveCompiledByte(Comp)
      end if
      
    case (cExp)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=exp(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cFloor)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=floor(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cIf)
      ! No optimization possible

    case (cLog)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) <= 0._DP) then
          call output_line('Invalid argument for LOG!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=log(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cLog10)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) <= 0._DP) then
          call output_line('Invalid argument for LOG!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=log10(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cMax)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=max(Comp%Immed(Comp%ImmedSize),Comp%Immed(Comp%ImmedSize-1))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cMin)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=min(Comp%Immed(Comp%ImmedSize),Comp%Immed(Comp%ImmedSize-1))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cSec)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        daux=cos(Comp%Immed(Comp%ImmedSize))
        if (daux == 0._DP) then
          call output_line('Invalid argument for SEC!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=1/daux
        call RemoveCompiledByte(Comp)
      end if

    case (cSin)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=sin(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cSinh)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=sinh(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if
      
    case (cSqrt)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) < 0._DP) then
          call output_line('Invalid argument for SQRT!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize)=sqrt(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cTan)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=tan(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cTanh)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=tanh(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

      !------------------------------------------------------------
      ! Misc
      !------------------------------------------------------------
    case (cImmed,cJump)
      ! No optimization needed

      !------------------------------------------------------------
      ! Operators
      !------------------------------------------------------------
    case (cNeg)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=-(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cAdd)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)+Comp%Immed(Comp%ImmedSize)
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cSub)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)-Comp%Immed(Comp%ImmedSize)
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cMul)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)*Comp%Immed(Comp%ImmedSize)
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cDiv)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) == 0._DP) then
          call output_line('Invalid argument for DIV!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)/Comp%Immed(Comp%ImmedSize)
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if
      
    case (cMod)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        if (Comp%Immed(Comp%ImmedSize) == 0._DP) then
          call output_line('Invalid argument for MOD!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        Comp%Immed(Comp%ImmedSize-1)=mod(Comp%Immed(Comp%ImmedSize-1),Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if
      
    case (cPow)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)**Comp%Immed(Comp%ImmedSize)
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cEqual)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) == Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cNEqual)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) /= Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if
      
    case (cLess)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) < Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cLessOrEq)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) <= Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cGreater)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) > Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cGreaterOrEq)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) >= Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if
      
    case (cAnd)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(DbleToLogc(Comp%Immed(Comp%ImmedSize-1)) .and.&
            DbleToLogc(Comp%Immed(Comp%ImmedSize)))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cOr)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .and.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) then
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(DbleToLogc(Comp%Immed(Comp%ImmedSize-1)) .or.&
            DbleToLogc(Comp%Immed(Comp%ImmedSize)))
        call RemoveCompiledImmediate(Comp)
        call RemoveCompiledByte(Comp)
        call RemoveCompiledByte(Comp)
      end if

    case (cNot)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=LogcToDble(.not.DbleToLogc(Comp%Immed(Comp%ImmedSize)))
        call RemoveCompiledByte(Comp)
      end if

      !------------------------------------------------------------
      ! Degrees-radians conversion
      !------------------------------------------------------------
    case (cDeg)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=RadToDeg(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if

    case (cRad)
      if (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) then
        Comp%Immed(Comp%ImmedSize)=DegToRad(Comp%Immed(Comp%ImmedSize))
        call RemoveCompiledByte(Comp)
      end if      
    end select
  end subroutine AddCompiledByte

  ! *****************************************************************************

!<subroutine>

  subroutine RemoveCompiledByte (Comp)

!<description>
    ! Remove last compiled byte from bytecode
!</description>

!<inputoutput>
    type (t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>
!</subroutine>
   
    Comp%ByteCode(Comp%ByteCodeSize) = 0
    Comp%ByteCodeSize = Comp%ByteCodeSize - 1
  end subroutine RemoveCompiledByte

  ! *****************************************************************************

!<subroutine>

  subroutine AddImmediate (Comp, immediate)

!<description>
    ! Add immediate
!</description>

!<input>
    ! Value of byte to be added
    real(DP), intent(in) :: immediate
!</input>

!<inputoutput>
    type (t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>
!</subroutine>
       
    Comp%ImmedSize = Comp%ImmedSize + 1
    Comp%Immed(Comp%ImmedSize) = immediate
  end subroutine AddImmediate

  ! *****************************************************************************

!<subroutine>

  subroutine RemoveCompiledImmediate (Comp)

!<description>
    ! Remove last compiled immediate from immediate stack
!</description>

!<inputoutput>
    type (t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    Comp%Immed(Comp%ImmedSize) = 0
    Comp%ImmedSize = Comp%ImmedSize - 1
  end subroutine RemoveCompiledImmediate

  ! *****************************************************************************

!<subroutine>

  subroutine AddFunctionOpcode (Comp, opcode)

!<description>
    ! Add compiled byte to bytecode
!</description>

!<input>
    ! Value of opcode to be added
    integer(is), intent(in) :: opcode
!</input>

!<inputoutput>
    type (t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>
!</subroutine>
    
    if (Comp%useDegreeConversion) then
      select case(opCode)
      case(cCos,Ccosh,cCot,cCsc,cSec,cSin,cSinh,cTan,cTanh)
        call AddCompiledByte(Comp,cRad)
      end select
    end if
    
    call AddCompiledByte(Comp,opcode)
    
    if (Comp%useDegreeConversion) then
      select case(opCode)
      case(cAcos,cAcosh,cAsinh,cAtanh,cAsin,cAtan,cAtan2)
        call AddCompiledByte(Comp, cDeg)
      end select
    end if
  end subroutine AddFunctionOpcode

  ! *****************************************************************************

!<function>

  function RealNum (str, ibegin, inext, error) result (res)

!<description>
    ! Get real number from string
    ! Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
!</description>

!<input>
    ! String
    character (LEN=*), intent(in) :: str
!</input>

!<output>
    ! OPTIONAL: Start position of real number
    integer, optional, intent(out) :: ibegin
    
    ! OPTIONAL: 1st character after real number
    integer, optional, intent(out) :: inext

    ! OPTIONAL: Error flag
    logical, optional, intent(out) :: error
!</output>

!<result>
    ! Real number
    real(DP) :: res
!</result>
!</function>
 
    ! local variables
    integer :: ib,in,istat
    logical :: Bflag,               & ! .T. at begin of number in str
               InMan,               & ! .T. in mantissa of number
               Pflag,               & ! .T. after 1st '.' encountered
               Eflag,               & ! .T. at exponent identifier 'eEdD'
               InExp,               & ! .T. in exponent of number
               DInMan,              & ! .T. if at least 1 digit in mant.
               DInExp,              & ! .T. if at least 1 digit in exp.
               err                    ! Local error flag

    
    Bflag=.true.; InMan=.false.; Pflag=.false.; Eflag=.false.; InExp=.false.
    DInMan=.false.; DInExp=.false.
    ib   = 1
    in   = 1
    do while (in <= len_trim(str))
      select case (str(in:in))
      case (' ') ! Only leading blanks permitted
        ib = ib+1
        if (InMan .or. Eflag .or. InExp) exit
      case ('+','-') ! Permitted only
        if     (Bflag) then           
          InMan=.true.; Bflag=.false. ! - at beginning of mantissa
        elseif (Eflag) then               
          InExp=.true.; Eflag=.false. ! - at beginning of exponent
        else
          exit ! - otherwise CALL sys_halt()
        endif
      case ('0':'9') ! Mark
        if     (Bflag) then           
          InMan=.true.; Bflag=.false. ! - beginning of mantissa
        elseif (Eflag) then               
          InExp=.true.; Eflag=.false. ! - beginning of exponent
        endif
        if (InMan) DInMan=.true. ! Mantissa contains digit
        if (InExp) DInExp=.true. ! Exponent contains digit
      case ('.')
        if     (Bflag) then
          Pflag=.true. ! - mark 1st appearance of '.'
          InMan=.true.; Bflag=.false. !   mark beginning of mantissa
        elseif (InMan .and..not.Pflag) then
          Pflag=.true. ! - mark 1st appearance of '.'
        else
          exit ! - otherwise CALL sys_halt()
        end if
      case ('e','E','d','D') ! Permitted only
        if (InMan) then
          Eflag=.true.; InMan=.false. ! - following mantissa
        else
          exit ! - otherwise CALL sys_halt()
        endif
      case DEFAULT
        exit ! CALL sys_halt() at all other characters
      end select
      in = in+1
    end do
    err = (ib > in-1) .or. (.not.DInMan) .or.&
        &((Eflag.or.InExp).and..not.DInExp)
    if (err) then
      res = 0.0_DP
    else
      read(str(ib:in-1),*,IOSTAT=istat) res
      err = istat /= 0
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
    if (present(error))  error  = err
  end function RealNum

  ! *****************************************************************************

!<function>

  recursive function FunctionSize (Func) result (fsize)

!<description>
    ! Return the size of the total function including external expressions
!</description>
    
!<input>
    ! Function string
    character (LEN=*),  intent(in) :: Func
!</input>

!<result>
    ! Size of function string
    integer :: fsize
!</result>
!</function>
    
    ! local variables
    integer :: ind,n
    character(LEN=1) :: c

    ! Determine size of given expression
    fsize=len_trim(Func)

    ! "Parse" string for externally defined expressions
    do ind=1,fsize
      c = Func(ind:ind)
      if (c == '@') then
        n = ExpressionIndex (Func(ind:))
        if (n == 0) then
          call output_line('Invalid expression!',&
                           OU_CLASS_WARNING,OU_MODE_STD,'FunctionSize')
          call sys_halt()
        end if
        fsize = fsize+FunctionSize(ExpressionStrings(n))
      end if
    end do
  end function FunctionSize

  ! *****************************************************************************

!<function>

  recursive function CompileExpression(Comp, Func, ind, Var, stopAtComma) result(ind2)

!<description>
    ! Compiles ','
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var

    ! OPTIONAL: stop at comma
    logical, intent(IN), optional :: stopAtComma
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: FuncLength
    
    FuncLength=len_trim(Func)
    
    ind2 = CompileOr(Comp,Func, ind, Var)
    if(ind2 > FuncLength) return

    if (present(stopAtComma)) then
      if (stopAtComma) return
    end if
    
    do while (Func(ind2:ind2) == ',')
      ind2 = CompileOr(Comp, Func, ind2+1, Var)
      if (ind2 > FuncLength) return
    end do
  end function CompileExpression

  ! *****************************************************************************

!<function>

  recursive function CompileOr(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles '|'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: FuncLength

    FuncLength=len_trim(Func)

    ind2 = CompileAnd(Comp, Func, ind, Var)
    if (ind2 > FuncLength) return

    do while(Func(ind2:ind2) == '|')
      ind2 = CompileAnd(Comp, Func, ind2+1, Var)

      call AddCompiledByte(Comp, cOr)
      Comp%StackPtr = Comp%StackPtr-1
      if (ind2 > FuncLength) return
    end do
  end function CompileOr

  ! *****************************************************************************

!<function>

  recursive function CompileAnd(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles '&'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: FuncLength

    FuncLength=len_trim(Func)
    
    ind2 = CompileComparison(Comp, Func, ind, Var)
    if (ind2 > FuncLength) return
    
    do while(Func(ind2:ind2) == '&')
      ind2 = CompileComparison(Comp, Func, ind2+1, Var)

      call AddCompiledByte(Comp, cAnd)
      Comp%StackPtr = Comp%StackPtr-1
      if (ind2 > FuncLength) return
    end do
  end function CompileAnd

  ! *****************************************************************************

!<function>

  recursive function CompileComparison(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles '=', '<' and '>'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer(is) :: opSize
    integer :: FuncLength

    FuncLength=len_trim(Func)
    
    ind2 = CompileAddition(Comp, Func, ind, Var)
    if (ind2 > FuncLength) return
    
    c=Func(ind2:ind2)
    do while(c == '=' .or. c == '<' .or. c == '>' .or. c == '!')
      opSize = merge(2,1,Func(ind2+1:ind2+1) == '=')
      ind2 = CompileAddition(Comp, Func, ind2+opSize, Var)

      select case(c)
      case('=')
        call AddCompiledByte(Comp, cEqual)

      case('<')
        call AddCompiledByte(Comp, merge(cLess,cLessOrEq,opSize==1))
        
      case('>')
        call AddCompiledByte(Comp, merge(cGreater,cGreaterOrEq,opSize==1))
        
      case('!')
        call AddCompiledByte(Comp, cNEqual)
      end select
      Comp%StackPtr = Comp%StackPtr-1
      
      if (ind2 > FuncLength) return
      c=Func(ind2:ind2)
    end do
  end function CompileComparison

  ! *****************************************************************************

!<function>

  recursive function CompileAddition(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles '+' and '-'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer :: FuncLength

    FuncLength=len_trim(Func)
    
    ind2 = CompileMult(Comp, Func, ind, Var)
    if (ind2 > FuncLength) return
    
    c=Func(ind2:ind2)
    do while(c == '+' .or. c == '-')
      ind2 = CompileMult(Comp, Func, ind2+1, Var)

      call AddCompiledByte(Comp, merge(cAdd,cSub,c == '+'))
      Comp%StackPtr = Comp%StackPtr-1
      
      if (ind2 > FuncLength) return
      c=Func(ind2:ind2)
    end do
  end function CompileAddition

  ! *****************************************************************************

!<function>

  recursive function CompileMult(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles '*', '/' and '%'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer :: FuncLength

    FuncLength=len_trim(Func)
    
    ind2 = CompileUnaryMinus(Comp, Func, ind, Var)
    if (ind2 > FuncLength) return
    
    c=Func(ind2:ind2)
    do while(c == '*' .or. c == '/' .or. c == '%')
      ind2 = CompileUnaryMinus(Comp, Func, ind2+1, Var)

      select case(c)
      case('*')
        call AddCompiledByte(Comp, cMul)
        
      case('/')
        call AddCompiledByte(Comp, cDiv)
        
      case('%')
        call AddCompiledByte(Comp, cMod)
        
      end select
      Comp%StackPtr = Comp%StackPtr-1

      if (ind2 > FuncLength) return
      c=Func(ind2:ind2)
    end do
  end function CompileMult

  ! *****************************************************************************

!<function>

  recursive function CompileUnaryMinus(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles unary '-'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer :: FuncLength

    FuncLength=len_trim(Func)

    c = Func(ind:ind)
    if (c == '-' .or. c == '!') then
      ind2 = ind+1
      if (ind2 > FuncLength) return
      ind2 = CompilePow(Comp, Func, ind2, Var)

      ! If we are negating a constant, negate the constant itself
      if (c == '-' .and. Comp%ByteCode(Comp%ByteCodeSize) == cImmed) then
        Comp%Immed(Comp%ImmedSize) = -Comp%Immed(Comp%ImmedSize)
        
        ! If we are negating a negation, we can remove both
      elseif (c == '-' .and. Comp%ByteCode(Comp%ByteCodeSize) == cNeg) then
        call RemoveCompiledByte(Comp)
        
      else
        call AddCompiledByte(Comp, merge(cNeg,cNot,c == '-'))

      end if
      return
    end if
    
    ind2 = CompilePow(Comp, Func, ind, Var)
  end function CompileUnaryMinus

  ! *****************************************************************************

!<function>

  recursive function CompilePow(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles '^'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: FuncLength

    FuncLength=len_trim(Func)

    ind2 = CompileElement(Comp, Func, ind, Var)

    if (ind2 > FuncLength) return

    do while(Func(ind2:ind2) == '^')
      ind2 = CompileUnaryMinus(Comp, Func, ind2+1, Var)

      call AddCompiledByte(Comp, cPow)
      Comp%StackPtr = Comp%StackPtr-1
      if (ind2 > FuncLength) return
    end do
  end function CompilePow

  ! *****************************************************************************

!<function>

  recursive function CompileElement(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles element
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    real(DP) :: rnum
    integer(is) :: n
    integer :: ind1,ind0,ib,in,requiredParams
    logical :: err

    ind1=ind; c=Func(ind1:ind1)
    if (c == '(') then
      ind1 = CompileExpression(Comp, Func, ind1+1, Var, .false.)
      ind2 = ind1+1   ! Func(ind1:ind1) is ')'
      return
    end if

    ! Check for numbers
    if (scan(c,'0123456789,') > 0) then
      rnum = RealNum (Func(ind1:),ib,in,err)
      if (err) then
        call output_line('Invalid number format!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'CompileElement')
        call sys_halt()
      end if
      call AddImmediate(Comp, rnum)
      call AddCompiledByte(Comp, cImmed)
      call incStackPtr(Comp)
      ind2 = ind1+in-1
      return
      
    else
      ! Then string must be function, variable or constant
      
      ! Check for mathematical functions
      n = MathFunctionIndex(Func(ind1:))
      if (n > 0) then
        ind2 = ind1+len_trim(Funcs(n))

        ! Check for IF-THEN-ELSE
        if (n == cIf) then
          ind2 = CompileIf(Comp, Func, ind2+1, Var)
          ! IF-THEN-ELSE cannot be vectorized, note that!
          Comp%isVectorizable = .false.
          return        
        end if

        requiredParams=MathFunctionParameters(n)        
        ind2 = CompileFunctionParameters(Comp, Func, ind2+1, Var, requiredParams)
        call AddFunctionOpcode(Comp, n)
        return
      end if
      
      ! Check for predefined constant
      n = ConstantIndex(Func(ind1:))
      if (n > 0) then
        ind2 = ind1+len_trim(Consts(n))+1
        call AddImmediate(Comp, Constvals(n))
        call AddCompiledByte(Comp, cImmed)
        call incStackPtr(Comp)
        return
      end if

      ! Check for predefined expressions
      n = ExpressionIndex(Func(ind1:))
      if (n > 0) then
        ind2 = ind1+len_trim(Expressions(n))+1

        ! Recursively compile the given expression
        ind0 = CompileExpression(Comp,ExpressionStrings(n), 1, Var)
        
        ! Proceed with compilation of mathematical function Func afterwards
        return
      end if

      ! Check for variables
      n = VariableIndex(Func(ind1:), Var, ib, in)
      if (n > 0) n = VarBegin+n-1
      call AddCompiledByte(Comp, n)
      call incStackPtr(Comp)
      ind2 = ind1+in-1
      return
    end if

    call output_line('An unexpected error occured!',&
                     OU_CLASS_WARNING,OU_MODE_STD,'CompileElement')
    call sys_halt()
  end function CompileElement

  ! *****************************************************************************

!<function>

  recursive function CompileFunctionParameters(Comp, Func, ind, Var,&
      & requiredParams) result(ind2)

!<description>
    ! Compiles function parameters
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var

    ! Number of required parameters
    integer, intent(IN) :: requiredParams
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: curStackPtr

    ind2 = ind
    if (requiredParams > 0) then
      
      curStackPtr = Comp%StackPtr
      ind2 = CompileExpression(Comp, Func, ind, Var, .false.)

      if (Comp%StackPtr /= curStackPtr+requiredParams) then
        call output_line('Illegal number of parameters to function!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'CompileFunctionParameters')
        call sys_halt()
      end if

      Comp%StackPtr = Comp%StackPtr-(requiredParams-1)

    else
      
      call incStackPtr(Comp)
      
    end if
    ind2=ind2+1
  end function CompileFunctionParameters

  ! *****************************************************************************

!<function>

  recursive function CompileIf(Comp, Func, ind, Var) result(ind2)

!<description>
    ! Compiles if()
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(IN) :: Func

    ! Begin position substring
    integer, intent(IN) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(INOUT) :: Comp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: FuncLength,curByteCodeSize,curByteCodeSize2,curImmedSize2

    FuncLength=len_trim(Func)

    ind2 = CompileExpression(Comp, Func, ind, Var, .true.) ! Condition branch
    if (ind2 > FuncLength) return

    if (Func(ind2:ind2) /= ',') then
      call output_line('Illegal number of parameters to function!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'CompileIf')
      call sys_halt()
    end if
    call AddCompiledByte(Comp, cIf)
    curByteCodeSize = Comp%ByteCodeSize
    call AddCompiledByte(Comp, 0_is) ! Jump index will be set below
    call AddCompiledByte(Comp, 0_is) ! Immed jump index will be set below
    Comp%StackPtr = Comp%StackPtr-1

    ind2 = CompileExpression(Comp, Func, ind2+1, Var, .true.) ! Then branch
    if (ind2 > FuncLength) return

    if (Func(ind2:ind2) /= ',') then
      call output_line('Illegal number of parameters to function!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'CompileIf')
      call sys_halt()
    end if
    call AddCompiledByte(Comp, cJump)
    curByteCodeSize2 = Comp%ByteCodeSize
    curImmedSize2 = Comp%ImmedSize
    call AddCompiledByte(Comp, 0_is) ! Jump index will be set below
    call AddCompiledByte(Comp, 0_is) ! Immed jump index will be set below
    Comp%StackPtr = Comp%StackPtr-1

    ind2 = CompileExpression(Comp, Func, ind2+1, Var, .true.) ! Else branch
    if (ind2 > FuncLength) return

    if (Func(ind2:ind2) /= ')') then
      call output_line('Illegal number of parameters to function!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'CompileIf')
      call sys_halt()
    end if
    
    ! Set jump indices
    if (associated(Comp%ByteCode)) then
      Comp%ByteCode(curByteCodeSize+1)  = curByteCodeSize2+2
      Comp%ByteCode(curByteCodeSize+2)  = curImmedSize2+1
      Comp%ByteCode(curByteCodeSize2+1) = Comp%ByteCodeSize
      Comp%ByteCode(curByteCodeSize2+2) = Comp%ImmedSize+1
    end if
    
    ind2=ind2+1
  end function CompileIf

  ! *****************************************************************************

!<function>

  elemental function DbleTOLogc(d) result(l)

!<description>
    ! This function transforms a Double into a Logical
!</description>

!<input>
    ! Double variable
    real(DP), intent(IN) :: d
!</input>

!<result>
    ! Logical variable
    logical :: l
!</result>
!</function>
    
    l=merge(.true.,.false.,abs(1-d)<=1e-12)
  end function DbleTOLogc

  ! *****************************************************************************

!<function>
  
  elemental function LogcToDble(l) result(d)

!<description>
    ! This function transforms a Logical into a Double
!</description>

!<input>
    ! Logical variable
    logical, intent(IN) :: l
!</input>

!<result>
    ! Double variable
    real(DP) :: d
!</result>
!</function>
      
    d=merge(1._DP,0._DP,l)
  end function LogcToDble

  ! *****************************************************************************

!<function>
  
  elemental function DegToRad(d) result(r)

!<description>
    ! This function converts DEG to RAD
!</description>

!<input>
    ! DEG
    real(DP), intent(IN) :: d
!</input>

!<result>
    ! RAD
    real(DP) :: r
!</result>
!</function>
      
    r=d*(3.141592653589793115997963468544185161590576171875_DP/180._DP)
  end function DegToRad

  ! *****************************************************************************

!<function>

  elemental function RadToDeg(r) result(d)
    
!<description>
    ! This function converts RAD to DEG
!</description>

!<input>
    ! RAD
    real(DP), intent(IN) :: r
!</input>

!<result>
    ! DEG
    real(DP) :: d
!</result>
!</function>
    
    d=r*(180._DP/3.141592653589793115997963468544185161590576171875_DP)
  end function RadToDeg
end module fparser
