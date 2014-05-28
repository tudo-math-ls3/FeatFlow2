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
!# ------------------------------------------------------------------------ \\
!# Copyright notice \\
!# ------------------------------------------------------------------------ \\
!# (a) This module is based on the "Fortran 90 function parser V1.1"
!#     written by Roland Schmehl < Roland.Schmehl@mach.uni-karlsruhe.de >
!#     The original Fortran90 source code is available from:
!#     http://itsextern.its.uni-karlsruhe.de/~schmehl/functionparser.html
!#
!#     However the "Fortran 90 function parser V1.1" only recognises the
!#     (single argument) Fortran 90 intrinsic functions abs, exp, log10,
!#     log, sqrt, sinh, cosh, tanh, sin, cos, tan, asin, acos, atan.
!#
!#     The function parser concept is based on a C++ class library written
!#     by Warp < warp@iki.fi > available from:
!#     http://warp.povusers.org/FunctionParser/
!#
!# (b) This FParser module is an extension of the `Fortran 90 function parser
!#     V1.1` which implements most of the features available in the `function
!#     parser library for C++ V2.8` written by Warp. The optimiser included
!#     in the C++ library and the recursive evaluation of functions by means
!#     of eval(...) is not implemented in this version.
!#
!# ------------------------------------------------------------------------ \\
!# Basic usage \\
!# ------------------------------------------------------------------------ \\
!#
!# Step 0 - Module Import \\
!# ---------------------- \\
!# In all program units where you want to use the function parser procedures
!# and variables you must import the module by:
!#
!# <code>
!#  use fparser
!# </code>
!#
!# This command imports only 6 public names: fparser_create, fparser_release,
!# fparser_parseFunction, fparser_evalFunction, fparser_ErrorMsg and EvalErrType
!# which are explained in the following. The remainder of the
!# module is hidden to the calling program.
!#
!# Step 1 - Initialization \\
!# ----------------------- \\
!# The parser module has to be initialised for the simultaneous evaluation of
!# n functions by calling the module subroutine initp one time in your Fortran
!# code:
!#
!# <code>
!#  call fparser_create (Parser, n)
!# </code>
!#
!# This allocates i=1,...,n internal data structures used by the byte-compiler
!# and subsequently by the bytecode-interpreter in the bytecode object Comp.
!#
!# Step 2 - Function parsing \\
!# ------------------------- \\
!# The i-th function string FuncStr is parsed (checked and compiled) into the
!# i-th bytecode by calling the module subroutine parsef:
!#
!# <code>
!#  call fparser_parseFunction (Parser, i, FuncStr, Var)
!# </code>
!#
!# The variable names as they appear in the string FuncStr have to be passed
!# in the one-dimensional string array Var (zero size of Var is acceptable).
!# The number of variables is implicitly passed by the dimension of this array.
!# For some notes on the syntax of the function string see below.
!#
!# Step 3 - Function evaluation \\
!# ---------------------------- \\
!# The i-th function value is evaluated for a specific set of variable values
!# by calling the module function evalf:
!#
!# <code>
!#  a = fparser_evalFunction (Parser, i, Val)
!# </code>
!#
!# The variable values are passed in the one-dimensional array Val which must
!# have the same dimension as array Var.
!#
!# ------------------------------------------------------------------------ \\
!# Error handling \\
!# ------------------------------------------------------------------------ \\
!#
!# An error in the function parsing step leads to a detailed error message
!# (Type and position of error) and program termination.
!#
!# An error during function evaluation returns a function value of 0.0 and
!# sets the error flag EvalErrType (part of the t_fparser derived type) to
!# a value > 0 (EvalErrType = 0 indicates no error). An error message from the
!# bytecode-interpreter can be obtained by calling the character function
!# fparser_ErrorMsg (Parser) with the parser object as an argument.
!#
!# ------------------------------------------------------------------------ \\
!# Function string syntax \\
!# ------------------------------------------------------------------------ \\
!#
!# Although they have to be passed as array elements of the same declared
!# length (Fortran 90 restriction), the variable names can be of arbitrary
!# actual length for the parser. Parsing for variables is case sensitive.
!#
!# The syntax of the function string is similar to the Fortran convention.
!# Mathematical Operators recognised are +, -, *, /, %, ** or alternatively
!# ^, whereas symbols for brackets must be (), [] or {}. Note that the
!# parser does not check if, e.g. ( is closed by ) or ]. At the moment,
!# different brackets may be used only to improve readability of the function
!# string.
!#
!# Operations are evaluated in the correct order:
!#
!# <verb>
!#  ()             expressions in brackets first
!#  -A             unary minus (or plus)
!#  A**B A^B       exponentiation (A raised to the power B)
!#  A*B  A/B  A%B  multiplication, division and modulo
!#  A+B  A-B       addition and subtraction
!#  A=B  A!=B  A < B  A <= B  A > B  A >= B
!#                 comparison between A and B (result is either 0 or 1)
!#  A&B            result is 1 if int(A) and int(B) differ from 0, else 0.
!#  A|B            result is 1 if int(A) or int(B) differ from 0, else 0.
!# </verb>
!#
!# The function string can contain integer or real constants. To be recognised
!# as explicit constants these must conform to the format
!#
!# <verb>
!#  [+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
!# </verb>
!#
!# where nnn means any number of digits. The mantissa must contain at least
!# one digit before or following an optional decimal point. Valid exponent
!# identifiers are 'e', 'E', 'd' or 'D'. If they appear they must be followed
!# by a valid exponent!
!#
!# Note that the function parser is case insensitive.
!# The following mathematical functions are supported
!#
!# <verb>
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
!# sign(A)   : Sign of A.
!# sqrt(A)   : Square root of A. Returns the value whose square is A.
!# tan(A)    : Tangent of A. Returns the tangent of the angle A, where A
!#             is measured in radians.
!# tanh(A)   : Same as tan() but for hyperbolic tangent.
!# rrand(A,B) : Reproducable pseudo-random number; B'th random number with
!#              Random-Seed A.
!# </verb>
!#
!# The parser also supports a number of standard constants in the function
!# string. All constants start with an underscore '_'. The following constants
!# are defined by default:
!#
!# <verb>
!# _PI       : Gives the number $pi$.
!# _EXP      : Gives the number $e$
!# _INFTY    : Gives the maximum possible number in double precision, defined
!#             in fsystem by SYS_INFINITY_DP.
!# </verb>
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
!#     -> Initialise the sub-system for function parsers
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
!# 7.) fparser_parseFunction = fparser_parseFunctionByName /
!#                             fparser_parseFunctionByNumber
!#     -> Parse function string and compile it into bytecode
!#
!# 8.) fparser_evalFunction = fparser_evalFuncScalarByName /
!#                            fparser_evalFuncScalarByNumber /
!#                            fparser_evalFuncBlockByName /
!#                            fparser_evalFuncBlockByNumber
!#     -> Evaluate precompiled bytecode
!#
!# 9.) fparser_ErrorMsg
!#     -> Get error message from function parser
!#
!# 10.) fparser_PrintByteCode = fparser_PrintByteCodeByName /
!#                              fparser_PrintByteCodeByNumber
!#      -> Print the bytecode stack (very technical!)
!#
!# 11.) fparser_parseFileForKeyword
!#      -> Parse input file for keyword
!#
!# 12.) fparser_getFunctionNumber
!#      -> Return the internal number of the function
!#
!# 13.) fparser_initPerfConfig
!#      -> Initialises the global performance configuration
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
!# 29.) evalFunctionScalar
!#      -> Evaluate function for scalar data
!#
!# 30.) evalFunctionBlock
!#      -> Evaluate function for multi-component data
!#
!# </purpose>
!##############################################################################

module fparser

  !$ use omp_lib
  use fsystem
  use genoutput
  use io
  use perfconfig
  use storage
  use stackInt

  implicit none

  private
  public :: t_fparser
  public :: fparser_initPerfConfig
  public :: fparser_init
  public :: fparser_done
  public :: fparser_parseFileForKeyword
  public :: fparser_defineConstant
  public :: fparser_defineExpression
  public :: fparser_create
  public :: fparser_release
  public :: fparser_parseFunction
  public :: fparser_evalFunction
  public :: fparser_evalFuncBlockByName2
  public :: fparser_evalFuncBlockByNumber2
  public :: fparser_ErrorMsg
  public :: fparser_PrintByteCode
  public :: fparser_getFunctionNumber

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  interface fparser_evalFunction
    module procedure fparser_evalFuncScalarByName
    module procedure fparser_evalFuncScalarByNumber
    module procedure fparser_evalFuncBlockByName
    module procedure fparser_evalFuncBlockByNumber
  end interface

  interface fparser_parseFunction
    module procedure fparser_parseFunctionByName
    module procedure fparser_parseFunctionByNumber
  end interface

  interface fparser_printByteCode
    module procedure fparser_printByteCodeByName
    module procedure fparser_printByteCodeByNumber
  end interface fparser_printByteCode

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<constants>

!<constantblock description="Global constants for parser">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of items to handle simultaneously when evaluating functions
#ifndef FPAR_NITEMSIM
  integer, parameter, public :: FPAR_NITEMSIM       = 256
#endif

  ! Length of string
  integer, parameter, public :: FPAR_STRLEN         = 2048

  ! Maximum number of predefined/user-defined constants
  integer, parameter, public :: FPAR_MAXCONSTS      = 128

  ! Maximum number of predefined/user-defined expressions
  integer, parameter, public :: FPAR_MAXEXPRESSIONS = 128

  ! Length of constant name
  integer, parameter, public :: FPAR_CONSTLEN       = 32

  ! Length of expression name
  integer, parameter, public :: FPAR_EXPRLEN        = 32

  ! Length of function name
  integer, parameter, public :: FPAR_FUNCLEN        = 32

  ! Length of variable name
  integer, parameter, public :: FPAR_VARLEN         = 32

!</constantblock>


!<constantblock description="types for parser expressions">

  ! Constant
  integer, parameter, public :: FPAR_CONSTANT   = 1

  ! Expression
  integer, parameter, public :: FPAR_EXPRESSION = 2

  ! Functions (including symbolic variables)
  integer, parameter, public :: FPAR_FUNCTION   = 3

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
                            cRrand       = 24, &
                            cAtan2       = 25, & ! --> last dyadic operator: .OP.(A,B)
                            cAbs         = 26, & ! <-- monadic operator: .OP.(A)
                            cAnint       = 27, &
                            cAint        = 28, &
                            cExp         = 29, &
                            cLog10       = 30, &
                            cLog         = 31
  integer(is), parameter :: cSqrt        = 32, &
                            cSinh        = 33, &
                            cCosh        = 34, &
                            cTanh        = 35, &
                            cSin         = 36, &
                            cCos         = 37, &
                            cTan         = 38, &
                            cCot         = 39, &
                            cAsinh       = 40, &
                            cAsin        = 41
  integer(is), parameter :: cAcosh       = 42, &
                            cAcos        = 43, &
                            cAtanh       = 44, &
                            cAtan        = 45, &
                            cCeil        = 46, &
                            cFloor       = 47, &
                            cCsc         = 48, &
                            cSec         = 49, &
                            cSign        = 50, & ! --> last monadic operator: .OP.(A)
                            VarBegin     = 51

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

  character (LEN=5), dimension(cIf:cSign), parameter :: Funcs = (/ 'if   ', &
                                                                   'min  ', &
                                                                   'max  ', &
                                                                   'rrand', &
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
                                                                   'sec  ', &
                                                                   'sign '/)

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
      SYS_INFINITY_DP/)

!</constantblock>


!<constantblock description="predefined expression names for parser; an at-sign '@' is automatically added">

  character(LEN=FPAR_CONSTLEN), dimension(1) :: PredefinedExpressions = (/ 'null      ' /)

!</constantblock>


!<constantblock description="predefined expressions for parser">

  character(LEN=FPAR_CONSTLEN), dimension(1) :: PredefinedExpressionvals = (/ '0         ' /)

!</constantblock>

!</constants>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<publicvars>

  ! Global number of predefined/user-defined constants
  integer, save :: nconstants = 0

  ! Global number of predefined/user-defined expressions
  integer, save :: nexpressions = 0

  ! Global constant names for parser
  character(LEN=FPAR_CONSTLEN), dimension(FPAR_MAXCONSTS), save :: CconstantName  = '     '

  ! Global constant values for parser
  real(DP), dimension(FPAR_MAXCONSTS), save :: DconstantValue = 0

  ! Global expression name for parser
  character(LEN=FPAR_EXPRLEN), dimension(FPAR_MAXCONSTS), save :: CexpressionName = ''

  ! Global expression string for parser
  character(LEN=FPAR_STRLEN), dimension(FPAR_MAXEXPRESSIONS), save :: CexpressionString

!</publicvars>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<types>

!<typeblock>

  ! Type block for storing all information of the function parser
  type t_fparser
    private

    ! Array of function parser components.
    ! Each component is used to handle one function string at a time
    type(t_fparserComponent), dimension(:), pointer :: Rcomp => null()

    ! Array of function names corresponding to the individual components.
    character(LEN=FPAR_FUNCLEN), dimension(:), pointer :: ScompName => null()

    ! Number of parser components
    integer :: ncomp = 0

    ! Maximum number of components
    integer :: nncomp = 0

  end type t_fparser

!</typeblock>

!<typeblock>

  ! Type block for storing the bytecode of the function parser for one component
  type t_fparserComponent
    private

    ! Size of bytecode
    integer :: ibytecodeSize = 0

    ! Size of immediates
    integer :: iimmedSize = 0

    ! Stack size
    integer :: istackSize = 0

    ! Stack pointer
    integer :: istackPtr = 0

    ! Use degree conversion DEG <-> RAD for some functions
    logical :: buseDegreeConversion = .false.

    ! Is vectorizable
    logical :: bisVectorizable = .true.

    ! Bytecode
    integer(is), dimension(:), pointer :: IbyteCode => null()

    ! Immediates
    real(DP), dimension(:), pointer :: Dimmed => null()
  end type t_fparserComponent
!</typeblock>

!</types>

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: fparser_perfconfig

  !************************************************************************

contains

  !****************************************************************************

!<subroutine>

  subroutine fparser_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      fparser_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(fparser_perfconfig)
      fparser_perfconfig%NELEMSIM = FPAR_NITEMSIM
    end if

  end subroutine fparser_initPerfConfig

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_init()

!<description>
    ! Initialise function parser
!</description>

!</subroutine>

    ! local variables
    integer :: i

    ! Initialise predefined constants
    do i = lbound(PredefinedConsts, 1),&
           ubound(PredefinedConsts, 1)
      call fparser_defineConstant(PredefinedConsts(i),&
                                  PredefinedConstvals(i))
    end do

    ! Initialise predefined expressions
    do i = lbound(PredefinedExpressions, 1),&
           ubound(PredefinedExpressions, 1)
      call fparser_defineExpression(PredefinedExpressions(i),&
                                    PredefinedExpressionvals(i))
    end do

  end subroutine fparser_init

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_done()

!<description>
    ! Release function parser
!</description>
!</subroutine>

    ! Reset constants
    nconstants     = 0
    Cconstantname  = '     '
    DconstantValue = 0._DP

    ! Reset expressions
    nexpressions      = 0
    CexpressionName   = ''
    CexpressionString = ''

  end subroutine fparser_done

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_parseFileForKeyword(rfparser, sfilename, ckeyword, itype)

!<description>
    ! Parse the file for the given keyword and make it a constant, a
    ! predefined expression or a function depending on the variable itype
!</description>

!<input>
    ! Name of parameter file to be parsed
    character(LEN=*), intent(in) :: sfilename

    ! Name of keyword to parser for
    character(LEN=*), intent(in) :: ckeyword

    ! Type of keyword: FPAR_CONSTANT, FPAR_EXPRESSION, FPAR_FUNCTION
    integer, intent(in) :: itype
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparser), intent(inout) :: rfparser
!</inputoutput>
!</subroutine>

    ! local variables
    character(FPAR_VARLEN), dimension(:), allocatable :: Svariables
    character(SYS_STRLEN) :: skeyword
    character(FPAR_CONSTLEN) :: sconstName
    character(FPAR_EXPRLEN) :: sexpressionName
    character(FPAR_FUNCLEN) :: sfunctionName
    character(FPAR_STRLEN) :: sdata,svalue,svariable
    integer :: iunit,ios,ipos,jpos,kpos,ivar,idatalen,icomp
    real(DP) :: dvalue

    ! Try to open the file
    call io_openFileForReading(sfilename, iunit, .true.)

    ! Oops...
    if (iunit .eq. -1) then
      call output_line('Unable to open input file!',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_parseFileForKeyword')
      call sys_halt()
    end if

    ! Read through the complete input file and look for global
    ! definitions of constants and fixed expressions
    ios = 0
    readline: do while(ios .eq. 0)

      ! Read next line in file
      call io_readlinefromfile(iunit, sdata, idatalen, ios)

      ! Check for keyword defconst or defexp
      ipos = scan(sdata(1:idatalen), ":")
      if (ipos .eq. 0) cycle

      call sys_tolower(sdata(1:max(1, ipos-1)), skeyword)
      if (trim(adjustl(skeyword)) .eq. ckeyword) then

        ! We found a keyword that will be applied to the parser
        select case(itype)
        case (FPAR_CONSTANT)

          ! Split the line into name and value
          jpos = scan(sdata(1:idatalen), "=")
          sconstName = trim(adjustl(sdata(ipos+1:jpos-1)))
          svalue = trim(adjustl(sdata(jpos+1:)))

          read(svalue,*) dvalue
          call fparser_defineConstant(sconstName, dvalue)

        case (FPAR_EXPRESSION)
          ! Split the line into name and value
          jpos = scan(sdata(1:idatalen), "=")
          sexpressionName = trim(adjustl(sdata(ipos+1:jpos-1)))
          svalue = trim(adjustl(sdata(jpos+1:)))

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
                  OU_CLASS_ERROR, OU_MODE_STD,'fparser_parseFileForKeyword')
              call sys_halt()
            end if

            ! Append line
            svalue(ipos:) = trim(adjustl(sdata(1:idatalen)))
          end do

          call fparser_defineExpression(sexpressionName, svalue)

        case (FPAR_FUNCTION)

          ! Split the line into name, expression and symbolic variables
          jpos = scan(sdata(1:idatalen), "=")
          sfunctionName = trim(adjustl(sdata(ipos+1:jpos-1)))
          svalue(1:) = sdata(jpos+1:idatalen)

          ! Check if function name already exists
          do icomp = 1, rfparser%ncomp
            if (trim(adjustl(sfunctionName)) .eq.&
                trim(rfparser%ScompName(icomp))) cycle readline
          end do

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
                  OU_CLASS_ERROR, OU_MODE_STD,'fparser_parseFileForKeyword')
              call sys_halt()
            end if

            ! Append line
            svalue(ipos:) = trim(adjustl(sdata(1:idatalen)))
          end do

          ! Extract the symbolic variables
          jpos = scan(svalue(1:ipos), ";", .true.)
          svariable = svalue(jpos+1:ipos)

          kpos = scan(svariable, ","); ivar = 0
          do while (kpos .ne. 0)
            ivar = ivar+1
            svariable = trim(adjustl(svariable(kpos+1:len_trim(svariable))))
            kpos = scan(svariable, ",")
          end do

          ! Allocate temporal memory
          allocate(Svariables(ivar+1))

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

          ! Parse the function
          call fparser_parseFunction(rfparser, sfunctionName, svalue(1:jpos-1), Svariables)

          ! Deallocate(temporal memory
          deallocate(Svariables)


        case default
          call output_line('Invalid type of expression!',&
              OU_CLASS_ERROR, OU_MODE_STD,'fparser_parseFileForKeyword')
          call sys_halt()
        end select

      end if
    end do readline

    ! Close file
    close (iunit)

  end subroutine fparser_parseFileForKeyword

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_defineConstant(sname, dvalue)

!<description>
    ! Define a new constant for the function parser.
    ! This subroutine checks if the given constant is already defined.
!</description>

!<input>
    ! Name of the constant
    character(LEN=FPAR_CONSTLEN), intent(in) :: sname

    ! Value of the constant
    real(DP), intent(in) :: dvalue
!</input>
!</subroutine>

    ! local variables
    character(LEN=len(sname)) :: sstring
    integer :: iconst

    ! Check if there is enough space
    if (nconstants .lt. FPAR_MAXCONSTS) then
      nconstants = nconstants+1
    else
      call output_line('No space left for definition of constant!',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_defineConstant')
      call sys_halt()
    end if

    ! Prepare constant
    call sys_tolower(sname, sstring)

    ! Check if constant is already defined
    do iconst = 1, nconstants-1
      if (CconstantName(iconst) .eq. sstring) then
        ! If it is already defined, then it must not have a different value
        if(DconstantValue(iconst) .ne. dvalue) then
          call output_line('Constant is already defined with different value!',&
              OU_CLASS_ERROR, OU_MODE_STD,'fparser_defineConstant')
          call sys_halt()
        else
          nconstants = nconstants-1
          return
        end if
      end if
    end do

    ! Apply constant value and constant name
    CconstantName(nconstants)  = sstring
    DconstantValue(nconstants) = dvalue

  end subroutine fparser_defineConstant

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_defineExpression(sname, svalue)

!<description>
    ! Define a new expression for the function parser.
    ! This subroutine checks if the given expression is already defined.
!</description>

!<input>
    ! Name of the expression
    character(LEN=FPAR_EXPRLEN), intent(in) :: sname

    ! String of the expression
    character(LEN=*), intent(in) :: svalue
!</input>
!</subroutine>

    ! local variables
    character(LEN=len(sname)) :: sexpression
    character(LEN=len(svalue)) :: sstring
    integer :: iexpression

    ! Check if there is enough space
    if (nexpressions .lt. FPAR_MAXEXPRESSIONS) then
      nexpressions = nexpressions+1
    else
      call output_line('No space left for definition of expression!',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_defineExpression')
      call sys_halt()
    end if

    ! Prepare expression string
    call sys_tolower(sname, sexpression)
    call sys_tolower(svalue, sstring)

    ! Replace human readable function names by 1-Char. format
    call Replace ('**','^ ', sstring)
    call Replace ('[','(',   sstring)
    call Replace (']',')',   sstring)
    call Replace ('{','(',   sstring)
    call Replace ('}',')',   sstring)

    ! Condense function string
    call RemoveSpaces (sstring)

    ! Check if expressions is already defined
    do iexpression = 1, nexpressions-1
      if (CexpressionName(iexpression) .eq. sexpression) then
        ! If it is already defined, then it must not have a different value
        if(CexpressionString(iexpression) .ne. sstring) then
          call output_line('Expression is already defined with different string!',&
              OU_CLASS_ERROR, OU_MODE_STD,'fparser_defineExpression')
          call sys_halt()
        else
          nexpressions = nexpressions-1
          return
        end if
      end if
    end do

    ! Apply expressions string and expression name
    CexpressionName(nexpressions)   = sexpression
    CexpressionString(nexpressions) = sstring

  end subroutine fparser_defineExpression

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_create (rfparser, nncomp)

!<description>
    ! Initialise function parser for nncomp functions.
!</description>

!<input>
    ! Number of functions
    integer, intent(in) :: nncomp
!</input>

!<output>
    ! Function parser object
    type(t_fparser), intent(out) :: rfparser
!</output>
!</subroutine>

    ! Set number of components
    rfparser%nncomp = nncomp
    rfparser%ncomp  = 0

    ! Allocate arrays
    allocate (rfparser%Rcomp(nncomp))
    allocate (rfparser%ScompName(nncomp))

  end subroutine fparser_create

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_release (rfparser)

!<description>
    ! Release function parser and all of its coponents
!</description>

!<inputoutput>
    ! Function parser
    type(t_fparser), intent(inout) :: rfparser
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: icomp

    ! Check that pointer is associated and return otherwise
    if (.not.associated(rfparser%Rcomp)) return

    ! Loop over all components and deallocate arrays
    do icomp = 1, rfparser%nncomp
      if (associated(rfparser%Rcomp(icomp)%IbyteCode)) deallocate(rfparser%Rcomp(icomp)%IbyteCode)
      if (associated(rfparser%Rcomp(icomp)%Dimmed))    deallocate(rfparser%Rcomp(icomp)%Dimmed)
    end do

    ! Deallocate memory
    deallocate(rfparser%Rcomp)
    deallocate(rfparser%ScompName)

    ! Reset scalar data
    rfparser%nncomp = 0
    rfparser%ncomp  = 0

  end subroutine fparser_release

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_parseFunctionByName (rfparser, scompName, sfunctionString,&
                                          Svariables, buseDegrees)

!<description>
    ! Parse function string sfuncStr and compile it into bytecode
!</description>

!<input>
    ! Function identifier
    character (LEN=*), intent(in) :: scompName

    ! Function string
    character (LEN=*), intent(in) :: sfunctionString

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Svariables

    ! OPTIONAL
    logical, intent(in), optional :: buseDegrees
!</input>

!<inputoutput>
    ! Function parser
    type (t_fparser), intent(inout) :: rfparser
!</inputoutput>
!</subroutine>

    ! Check if there is space for a new component
    if (rfparser%ncomp .eq. rfparser%nncomp) then
      call output_line('No free components left!',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_parseFunctionByName')
      call sys_halt()
    end if

    ! Increase component counter
    rfparser%ncomp = rfparser%ncomp+1

    ! Store function name
    call sys_tolower(trim(adjustl(scompName)), rfparser%ScompName(rfparser%ncomp))

    ! Parse function
    call fparser_parseFunction (rfparser, rfparser%ncomp, sfunctionString,&
                                Svariables, buseDegrees)

  end subroutine fparser_parseFunctionByName

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_parseFunctionByNumber (rfparser, icomp, sfunctionString,&
                                            Svariables, buseDegrees)

!<description>
    ! Parse ith function string sfuncStr and compile it into bytecode
!</description>

!<input>
    ! Function identifier
    integer, intent(in) :: icomp

    ! Function string
    character (LEN=*), intent(in) :: sfunctionString

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Svariables

    ! OPTIONAL
    logical, intent(in), optional :: buseDegrees
!</input>

!<inputoutput>
    ! Function parser
    type (t_fparser), intent(inout) :: rfparser
!</inputoutput>
!</subroutine>

    ! local variables
    character (LEN=len(sfunctionString)) :: sstring

    ! Check if component is valid
    if (icomp .lt. 1 .or. icomp .gt. rfparser%nncomp) then
      call output_line('Component number is out of range',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_parseFunctionByNumber')
      call sys_halt()
    end if

    ! Local copy of function string
    sstring = sfunctionString

    ! Replace human readable function names by 1-Char. format
    call Replace ('**','^ ', sstring)
    call Replace ('[','(',   sstring)
    call Replace (']',')',   sstring)
    call Replace ('{','(',   sstring)
    call Replace ('}',')',   sstring)

    ! Condense function string
    call RemoveSpaces (sstring)

    ! Check for valid syntax; this prevents the bytecode compiler
    ! from running into endless loops or other problems
    call CheckSyntax (sstring, Svariables)

    ! Check if conversion to degrees is required
    if (present(buseDegrees))&
        rfparser%Rcomp(icomp)%buseDegreeConversion = buseDegrees

    ! If syntax is correct, then compile into bytecode
    call Compile (rfparser%Rcomp(icomp), sstring, Svariables)

  end subroutine fparser_parseFunctionByNumber

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_evalFuncScalarByName (rfparser, scompName, Dvalue, dresult)

!<description>
    ! Evaluate bytecode of function named scompName for the values
    ! passed in array Cval(:). Note that this function is a wrapper for
    ! the working routine evalFunctionScalar. It is used to adjust the
    ! dimensions of the global stack memory if required.
!</description>

!<input>
    ! Function parser
    type (t_fparser), intent(in) :: rfparser

    ! Function name
    character(LEN=*), intent(in) :: scompName

    ! Variable values
    real(DP), dimension(:), intent(in) :: Dvalue
!</input>

!<output>
    ! Evaluated function
    real(DP), intent(out)  :: dresult
!</output>
!</subroutine>

    ! local variables
    integer :: icomp

    ! Lookup function by name
    icomp = fparser_getFunctionNumber(rfparser, scompName)

    ! Evaluate function by number
    call fparser_evalFunction (rfparser, icomp, Dvalue, dresult)

  end subroutine fparser_evalFuncScalarByName

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_evalFuncBlockByName (rfparser, scompName, idim, DValueBlock,&
                                          Dresult, DvalueScalar, rperfconfig)

!<description>
    ! Evaluate bytecode of component icomp for the array of values passed in
    ! array DvalueBlock(:,:). Note that this function is a wrapper for the working
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
    ! and DvalueScalar should be used, then the first variables must stored as
    ! blocks whereas the last variables can be scalar. This sound slightly
    ! complicated but here is an example:
    !
    ! Suppose you want to evaluate a function f=f(x,y,t). You know that x,y
    ! corresponds to the coordinate vector and t denotes the time. Then
    ! you should order your variables according to [x,y,t]. If the function
    ! should be evaluated for a set of variables then DvalueBlock=[x,y] and
    ! DvalueScalar=[t] works fine.
!</description>

!<input>
    ! Function parser
    type (t_fparser), intent(in) :: rfparser

    ! Function name
    character(LEN=*), intent(in) :: scompName

    ! Orientation of the stored values
    ! idim =1 : DvalueBlock is organised as (x1:xN),(y1:yN),...
    ! idim =2 : DvalueBlock is organised as (x1,y1),(x2,y2),...,(xN,yN)
    integer, intent(in) :: idim

    ! Variable values (must have the same dimension as Dresult)
    real(DP), dimension(:,:), intent(in) :: DvalueBlock

    ! Variable values. This is a vector of scalar variables
    ! which is the same for all components of Res, e.g. the time variable.
    real(DP), dimension(:), intent(in), optional :: DvalueScalar

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
    ! Evaluated function
    real(DP), dimension(:), intent(out) :: Dresult
!</output>
!</subroutine>

    ! local variables
    integer :: icomp

    ! Lookup function by name
    icomp = fparser_getFunctionNumber(rfparser, scompName)

    ! Evaluate function by number
    call fparser_evalFunction (rfparser, icomp, idim, DvalueBlock,&
                               Dresult, DvalueScalar, rperfconfig)

  end subroutine fparser_evalFuncBlockByName

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_evalFuncBlockByName2 (rfparser, scompName, n1, n2, DValueBlock,&
                                           n3, Dresult, DvalueScalar, rperfconfig)

!<description>
    ! Evaluate bytecode of component icomp for the array of values passed in
    ! array DvalueBlock(:,:). Note that this function is a wrapper for the working
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
    ! and DvalueScalar should be used, then the first variables must stored as
    ! blocks whereas the last variables can be scalar. This sound slightly
    ! complicated but here is an example:
    !
    ! Suppose you want to evaluate a function f=f(x,y,t). You know that x,y
    ! corresponds to the coordinate vector and t denotes the time. Then
    ! you should order your variables according to [x,y,t]. If the function
    ! should be evaluated for a set of variables then DvalueBlock=[x,y] and
    ! DvalueScalar=[t] works fine.
!</description>

!<input>
    ! Function parser
    type (t_fparser), intent(in) :: rfparser

    ! Function name
    character(LEN=*), intent(in) :: scompName

    ! Array dimensions
    integer, intent(in) :: n1,n2,n3

    ! Variable values (must have the same dimension as Dresult)
    real(DP), dimension(n1,n2), intent(in) :: DvalueBlock

    ! Variable values. This is a vector of scalar variables
    ! which is the same for all components of Res, e.g. the time variable.
    real(DP), dimension(:), intent(in), optional :: DvalueScalar

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
    ! Evaluated function
    real(DP), dimension(n3), intent(out) :: Dresult
!</output>
!</subroutine>

    ! local variables
    integer :: icomp

    ! Lookup function by name
    icomp = fparser_getFunctionNumber(rfparser, scompName)

    ! Evaluate function by number
    if (n1 .eq. n3) then
      call fparser_evalFunction (rfparser, icomp, 1, DvalueBlock,&
                                 Dresult, DvalueScalar, rperfconfig)
    elseif (n2 .eq. n3) then
      call fparser_evalFunction (rfparser, icomp, 2, DvalueBlock,&
                                 Dresult, DvalueScalar, rperfconfig)
    else
      call output_line('Invalid array dimensions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fparser_evalFuncBlockByName2')
      call sys_halt()
    end if

  end subroutine fparser_evalFuncBlockByName2

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_evalFuncScalarByNumber (rfparser, icomp, Dvalue, dresult)

!<description>
    ! Evaluate bytecode of component icomp for the values passed in array
    ! Dvalue(:). Note that this function is a wrapper for the working routine
    ! evalFunctionScalar. It is used to adjust the dimensions of the global
    ! stack memory if required.
!</description>

!<input>
    ! Function parser
    type (t_fparser), intent(in) :: rfparser

    ! Function identifier
    integer, intent(in) :: icomp

    ! Variable values
    real(DP), dimension(:), intent(in) :: Dvalue
!</input>

!<output>
    ! Evaluated function
    real(DP), intent(out)  :: dresult
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:), allocatable :: Dstack
    integer :: EvalErrType

    ! Check if component is valid
    if (icomp .lt. 1 .or. icomp .gt. rfparser%nncomp) then
      call output_line('Component number is out of range',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_evalFuncScalarByNumber')
      call sys_halt()
    end if

    ! Allocate temporal memory
    allocate(Dstack(rfparser%Rcomp(icomp)%istackSize+1))

    ! Invoke working routine
    call evalFunctionScalar(rfparser%Rcomp(icomp), Dstack,&
                            Dvalue, EvalErrType, dresult)

    ! Deallocate temporal memory
    deallocate(Dstack)

    ! Check if evaluation was successful
    if (EvalErrType .ne. 0) then
      call output_line('An error occured during function evaluation!',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_evalFuncScalarByNumber')
      call sys_halt()
    end if

  end subroutine fparser_evalFuncScalarByNumber

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_evalFuncBlockByNumber (rfparser, icomp, idim, DvalueBlock,&
                                            Dresult, DvalueScalar, rperfconfig)

!<description>
    ! Evaluate bytecode of component icomp for the array of values passed in
    ! array DvaluelBlock(:,:). Note that this function is a wrapper for the working
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
    ! Suppose you want to evaluate a function f=f(x,y,t). You know that x,y
    ! corresponds to the coordinate vector and t denotes the time. Then
    ! you should order your variables according to [x,y,t]. If the function
    ! should be evaluated for a set of variables then DvalueBlock=[x,y] and
    ! DvalueScalar=[t] works fine.
!</description>

!<input>
    ! Function parser
    type (t_fparser), intent(in) :: rfparser

    ! Function identifier
    integer, intent(in) :: icomp

    ! Orientation of the stored values
    ! idim =1 : DvalueBlock is organised as (x1:xN),(y1:yN),...
    ! idim =2 : DvalueBlock is organised as (x1,y1),(x2,y2),...,(xN,yN)
    integer, intent(in) :: idim

    ! Variable values (must have the same dimension as Dresult)
    real(DP), dimension(:,:), intent(in) :: DvalueBlock

    ! Variable values. This is a vector of scalar variables
    ! which is the same for all components of Res, e.g. the time variable.
    real(DP), dimension(:), intent(in), optional :: DvalueScalar

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
    ! Evaluated function
    real(DP), dimension(:), intent(out) :: Dresult
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), allocatable :: Dstack
    real(DP), dimension(:), allocatable :: DvalueTemp
    integer :: iValSet,iValMax,nvalue,iblockSize,isizeValueScalar,isizeValueBlock
    integer :: EvalErrType

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => fparser_perfconfig
    end if

    ! Check if component is valid
    if (icomp .lt. 1 .or. icomp .gt. rfparser%nncomp) then
      call output_line('Component number is out of range',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_evalFuncBlockByNumber')
      call sys_halt()
    end if

    ! Get total number of variable sets
    nvalue = size(DvalueBlock, idim)

    ! Initialise error flag
    EvalErrType = 0

    ! Check if the compiled function is vectorizable
    if (rfparser%Rcomp(icomp)%bisVectorizable) then
      ! ...ok, vectorization of the bytecode is admissible.

      ! What is the organization of ValBlock(:,:)
      if (idim .eq. 1) then
        !$omp parallel default(shared)&
        !$omp private(Dstack,iValMax,iblockSize)&
        !$omp reduction(max:EvalErrType)

        ! Allocate temporal memory
        allocate(Dstack(p_rperfconfig%NITEMSIM,rfparser%Rcomp(icomp)%iStackSize+1))

        !$omp do schedule(static,1)
        do iValSet = 1, nvalue, p_rperfconfig%NITEMSIM

          ! Initialization
          iValMax    = min(iValSet+p_rperfconfig%NITEMSIM-1, nvalue)
          iblockSize = iValMax-iValSet+1

          ! Invoke working routine
          call evalFunctionBlock(rfparser%Rcomp(icomp), iblockSize, Dstack,&
                                 DvalueBlock(iValSet:iValMax,:), idim, EvalErrType,&
                                 Dresult(iValSet:iValMax), DvalueScalar)
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(Dstack)
        !$omp end parallel
      else
        !$omp parallel default(shared)&
        !$omp private(Dstack,iValMax,iblockSize)&
        !$omp reduction(max:EvalErrType)

        ! Allocate temporal memory
        allocate(Dstack(p_rperfconfig%NITEMSIM,rfparser%Rcomp(icomp)%iStackSize+1))

        !$omp do schedule(static,1)
        do iValSet = 1, nvalue, p_rperfconfig%NITEMSIM

          ! Initialization
          iValMax    = min(iValSet+p_rperfconfig%NITEMSIM-1, nvalue)
          iblockSize = iValMax-iValSet+1

          ! Invoke working routine
          call evalFunctionBlock(rfparser%Rcomp(icomp), iblockSize, Dstack,&
                                 DvalueBlock(:, iValSet:iValMax), idim, EvalErrType,&
                                 Dresult(iValSet:iValMax), DvalueScalar)
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(Dstack)
        !$omp end parallel
      end if

    else   ! The compiled function cannot be vectorised

      ! Allocate temporal memory
      allocate(Dstack(rfparser%Rcomp(icomp)%iStackSize+1,1))

      ! The compiled bytecode cannot be vectorised. Hence, evaluate the function
      ! separately for each set of variables. Here, the organization of the array
      ! DvalBlock(:,:) is important. Moreover, if the optional parameter DvalScalar is
      ! given, then we have to combine those variables from DvalBlock and DvalScalar.

      if (present(DvalueScalar)) then

        ! Allocate auxiliary array
        isizeValueBlock  = size(DvalueBlock,3-idim)
        isizeValueScalar = size(DvalueScalar)
        allocate(DvalueTemp(isizeValueBlock+isizeValueScalar))

        if (idim .eq. 1) then
          do iValSet = 1, nvalue

            DvalueTemp(1:isizeValueBlock)  = DvalueBlock(iValSet,1:isizeValueBlock)
            DvalueTemp(isizeValueBlock+1:) = DvalueScalar

            ! Invoke working routine
            call evalFunctionScalar(rfparser%Rcomp(icomp), Dstack(:,1),&
                                    DvalueTemp, EvalErrType, Dresult(iValSet))
          end do
        else
          do iValSet = 1, nvalue

            DvalueTemp(1:isizeValueBlock)  = DvalueBlock(:,iValSet)
            DvalueTemp(isizeValueBlock+1:) = DvalueScalar

            ! Invoke working routine
            call evalFunctionScalar(rfparser%Rcomp(icomp), Dstack(:,1),&
                                    DvalueTemp, EvalErrType, Dresult(iValSet))
          end do
        end if

        ! Deallocate auxiliary array
        deallocate(DvalueTemp)

      else

        if (idim .eq. 1) then
          do iValSet = 1, nvalue

            ! Invoke working routine
            call evalFunctionScalar(rfparser%Rcomp(icomp), Dstack(:,1),&
                                    DvalueBlock(iValSet,:), EvalErrType,&
                                    Dresult(iValSet))
          end do
        else
          do iValSet = 1, nvalue

            ! Invoke working routine
            call evalFunctionScalar(rfparser%Rcomp(icomp), Dstack(:,1),&
                                    DvalueBlock(:, iValSet), EvalErrType,&
                                    Dresult(iValSet))
          end do
        end if

      end if

      ! Deallocate temporal memory
      deallocate(Dstack)

    end if

    ! Check if evaluation was successful
    if (EvalErrType .ne. 0) then
      call output_line('An error occured during function evaluation!',&
          OU_CLASS_ERROR, OU_MODE_STD,'fparser_evalFuncBlockByNumber')
      call sys_halt()
    end if

  end subroutine fparser_evalFuncBlockByNumber

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_evalFuncBlockByNumber2 (rfparser, icomp, n1, n2, DvalueBlock,&
                                             n3, Dresult, DvalueScalar, rperfconfig)

!<description>
    ! Evaluate bytecode of component icomp for the array of values passed in
    ! array DvaluelBlock(:,:). Note that this function is a wrapper for the working
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
    ! Suppose you want to evaluate a function f=f(x,y,t). You know that x,y
    ! corresponds to the coordinate vector and t denotes the time. Then
    ! you should order your variables according to [x,y,t]. If the function
    ! should be evaluated for a set of variables then DvalueBlock=[x,y] and
    ! DvalueScalar=[t] works fine.
!</description>

!<input>
    ! Function parser
    type (t_fparser), intent(in) :: rfparser

    ! Function identifier
    integer, intent(in) :: icomp

    ! Array dimensions
    integer, intent(in) :: n1,n2,n3

    ! Variable values (must have the same dimension as Dresult)
    real(DP), dimension(n1,n2), intent(in) :: DvalueBlock

    ! Variable values. This is a vector of scalar variables
    ! which is the same for all components of Res, e.g. the time variable.
    real(DP), dimension(:), intent(in), optional :: DvalueScalar

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
    ! Evaluated function
    real(DP), dimension(n3), intent(out) :: Dresult
!</output>
!</subroutine>

    ! Evaluate function by number
    if (n1 .eq. n3) then
      call fparser_evalFunction (rfparser, icomp, 1, DvalueBlock,&
                                 Dresult, DvalueScalar, rperfconfig)
    elseif (n2 .eq. n3) then
      call fparser_evalFunction (rfparser, icomp, 2, DvalueBlock,&
                                 Dresult, DvalueScalar, rperfconfig)
    else
      call output_line('Invalid array dimensions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fparser_evalFuncBlockByNumber2')
      call sys_halt()
    end if

  end subroutine fparser_evalFuncBlockByNumber2

  ! *****************************************************************************

!<function>

  function fparser_ErrorMsg (EvalErrType) result (smessage)

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
    integer, intent(in) :: EvalErrType
!</input>

!<result>
    ! Error messages
    character (LEN=len(m)) :: smessage
!</result>
!</function>

    if (EvalErrType .lt. 1 .or. EvalErrType .gt. size(m)) then
      smessage = ''
    else
      smessage = m(EvalErrType)
    endif

  end function fparser_ErrorMsg

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_PrintByteCodeByName(rfparser, scompName)

!<description>
    ! Print the compiled bytecode stack
!</description>

!<input>
    ! Function parser
    type(t_fparser), intent(in) :: rfparser

    ! Function name
    character(LEN=*), intent(in) :: scompName
!</input>
!</subroutine>

    ! local variables
    integer :: icomp

    ! Lookup function by name
    icomp = fparser_getFunctionNumber(rfparser, scompName)

    ! Print bytecode
    call fparser_PrintByteCodeByNumber(rfparser, icomp)

  end subroutine fparser_PrintByteCodeByName

  ! *****************************************************************************

!<subroutine>

  subroutine fparser_PrintByteCodeByNumber(rfparser, icomp)

!<description>
    ! Print the compiled bytecode stack
!</description>

!<input>
    ! Function parser
    type(t_fparser), intent(in) :: rfparser

    ! Function identifier
    integer, intent(in) :: icomp
!</input>
!</subroutine>

    ! local variables
    type(t_fparserComponent), pointer :: p_Comp
    character(LEN=5) :: n
    integer :: iinstPtr, idataPtr, istackPtr, nparams
    integer(is) :: iopCode

    nparams   = 1
    idataPtr  = 1
    istackPtr = 0
    iinstPtr  = 0
    p_Comp    => rfparser%Rcomp(icomp)

    do while(iinstPtr .lt. p_Comp%ibytecodeSize)
      iinstPtr = iinstPtr+1

      write(*,FMT='(I8.8,1X,":",1X)', ADVANCE="NO") iinstPtr

      iopCode = p_Comp%IbyteCode(iinstPtr)
      select case(iopCode)

      case(cIf)
        write(*,FMT='(A,1X,T10,I8.8)') "jz", p_Comp%IbyteCode(iinstPtr+1)+1
        iinstPtr = iinstPtr+2

      case(cJump)
        write(*,FMT='(A,1X,T10,I8.8)') "jump", p_Comp%IbyteCode(iinstPtr+1)+1
        iinstPtr = iinstPtr+2

      case(cImmed)
        write(*,FMT='(A,1X,T10,G16.8)') "push", p_Comp%Dimmed(idataPtr)
        idataPtr = idataPtr+1

      case default
        if (iopCode .lt. VarBegin) then
          select case(iopCode)
          case(cNEG);         n = "neg"
          case(cADD);         n = "add"
          case(cSUB);         n = "sub"
          case(cMUL);         n = "mul"
          case(cDIV);         n = "div"
          case(cMOD);         n = "mod"
          case(cPOW);         n = "pow"
          case(cEqual);       n = "eq"
          case(cNEqual);      n = "ne"
          case(cLess);        n = "lt"
          case(cLessOrEq);    n = "le"
          case(cGreater);     n = "gt"
          case(cGreaterOrEq); n = "ge"
          case(cAND);         n = "and"
          case (cOR);         n = "or"
          case(cNOT);         n = "not"
          case(cDEG);         n = "deg"
          case(cRAD);         n = "rad"

          case default
            n       = Funcs(iopCode)
            nparams = MathFunctionParameters(iopCode)
          end select
          write(*,FMT='(A,T10,A,"  (",I1,") ")') trim(n), "Par", nparams

        else
          write(*,FMT='(A,T10,A,1X,I4.4)') "push", "Var", iopCode-VarBegin+1
        end if

      end select
    end do

  end subroutine fparser_PrintByteCodeByNumber

  ! *****************************************************************************

!<function>

  function fparser_getFunctionNumber (rfparser, scompName, bquiet) result(icomp)

!<description>
    ! This function returns the internal number of the component which
    ! correspones to the function with name scompName
!</description>

!<input>
    ! Function parser
    type(t_fparser), intent(in) :: rfparser

    ! Function name
    character(len=*), intent(in) :: scompName

    ! OPTIONAL: Specifies whether a warning should be printed when released an
    ! empty vector (bquiet = .false.) or whether to remain silent in this case.
    ! If not specified, bquiet = .false. is used.
    logical, optional, intent(in) :: bquiet
!</input>

!<result>
    ! Function identifier
    integer :: icomp
!</result>
!</function>

    ! local variable
    logical :: bwarn
    character(len=len(scompName)) :: sname

    ! Convert to lower case
    call sys_tolower(scompName, sname)

    ! Lookup function
    do icomp = 1, rfparser%ncomp
      if (trim(adjustl(sname)) .eq. trim(rfparser%ScompName(icomp))) return
    end do

    ! If we end up here, then the function is not available
    icomp = 0

    ! Shout or shut up?
    bwarn = .true.
    if(present(bquiet)) bwarn = .not. bquiet

    if (bwarn) call output_line('Function is not available',&
        OU_CLASS_WARNING, OU_MODE_STD,'fparser_getFunctionNumber')

  end function fparser_getFunctionNumber

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<subroutine>

  recursive subroutine CheckSyntax (sfunctionString, Svariables)

!<description>
    ! Check syntax of function string, returns 0 if syntax is ok
!</description>

!<input>
    ! Function string without spaces
    character (LEN=*), intent(in) :: sfunctionString

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Svariables
!</input>
!</subroutine>

    ! local avriables
    type(t_stackInt) :: rstack
    integer(is) :: n,iopSize
    character (LEN=1) :: c
    real(DP) :: dnumber
    logical :: berror
    integer :: ifunctionIndex,ifunctionIndex2
    integer :: iparenthCount,ib,in,ifunctionLength,idummy

    ! Initialization
    ifunctionIndex  = 1
    iparenthCount   = 0
    ifunctionLength = len_trim(sfunctionString)
    call stack_create(rstack, max(5, int(ifunctionLength/4._DP)))

    do
      if (ifunctionIndex .gt. ifunctionLength) then
        call output_line('Invalid function string '//&
              trim(adjustl(sfunctionString))//' !',&
            OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
        call sys_halt()
      end if
      c = sfunctionString(ifunctionIndex:ifunctionIndex)

      ! Check for valid operand (must appear)

      ! Check for leading - or !
      if (c .eq. '-' .or. c .eq. '!') then
        ifunctionIndex = ifunctionIndex+1
        if (ifunctionIndex .gt. ifunctionLength) then
          call output_line('Premature end of string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        c = sfunctionString(ifunctionIndex:ifunctionIndex)
      end if

      ! Check for math function
      n = MathFunctionIndex (sfunctionString(ifunctionIndex:))
      if (n .gt. 0) then
        ! Math function found
        ifunctionIndex = ifunctionIndex+len_trim(Funcs(n))
        if (ifunctionIndex > ifunctionLength) then
          call output_line('Premature end of string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        c = sfunctionString(ifunctionIndex:ifunctionIndex)
        if (c .ne. '(') then
          call output_line('Expecting ( after function '//sfunctionString(ifunctionIndex:)//'!',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        ifunctionIndex2 = ifunctionIndex+1
        if (ifunctionIndex2 .gt. ifunctionLength) then
          call output_line('Premature end of string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        if (sfunctionString(ifunctionIndex2:ifunctionIndex2) .eq. ')') then
          ifunctionIndex = ifunctionIndex2+1
          if (ifunctionIndex .gt. ifunctionLength) then
            call output_line('Premature end of string '//&
              trim(adjustl(sfunctionString))//' !',&
                OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
            call sys_halt()
          end if
          c = sfunctionString(ifunctionIndex:ifunctionIndex)

          ! Ugly, but other methods would just be uglier ...
          goto 999
        end if

        ! Push counter for parenthesss to stack
        call stack_push(rstack, iparenthCount+1)
      end if

      ! Check for opening parenthesis
      if (c .eq. '(') then
        iparenthCount = iparenthCount+1
        ifunctionIndex = ifunctionIndex+1
        if (ifunctionIndex .gt. ifunctionLength) then
          call output_line('Premature end of string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        if (sfunctionString(ifunctionIndex:ifunctionIndex) .eq. ')') then
          call output_line('Empty parantheses in string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        cycle
      end if

      ! Check for number
      if (scan(c,'0123456789.') .gt. 0) then
        dnumber = RealNum (sfunctionString(ifunctionIndex:), ib, in, berror)
        if (berror) then
          call output_line('Invalid number format in string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        ifunctionIndex = ifunctionIndex+in-1
        if (ifunctionIndex > ifunctionLength) exit
        c = sfunctionString(ifunctionIndex:ifunctionIndex)
      elseif (c .eq. '_') then
        ! Check for constant
        n = ConstantIndex (sfunctionString(ifunctionIndex:))
        if (n .eq. 0) then
          call output_line('Invalid constant in string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        ifunctionIndex = ifunctionIndex+len_trim(CconstantName(n))+1
        if (ifunctionIndex > ifunctionLength) exit
        c = sfunctionString(ifunctionIndex:ifunctionIndex)
      elseif (c .eq. '@') then
        ! Check for expression
        n = ExpressionIndex (sfunctionString(ifunctionIndex:))
        if (n .eq. 0) then
          call output_line('Invalid expression in string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        call CheckSyntax(CexpressionString(n), Svariables)
        ifunctionIndex = ifunctionIndex+len_trim(CexpressionName(n))+1
        if (ifunctionIndex .gt. ifunctionLength) exit
        c = sfunctionString(ifunctionIndex:ifunctionIndex)
      else
        ! Check for variable
        n = VariableIndex (sfunctionString(ifunctionIndex:), Svariables, ib, in)
        if (n .eq. 0) then
          call output_line('Invalid element in string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        ifunctionIndex = ifunctionIndex+in-1
        if (ifunctionIndex .gt. ifunctionLength) exit
        c = sfunctionString(ifunctionIndex:ifunctionIndex)
      end if

      ! Check for closing parenthesis
      do while (c .eq. ')')
        if (.not.stack_empty(rstack)) then
          call stack_top(rstack, idummy)
          if(idummy .eq. iparenthCount) call stack_pop(rstack, idummy)
        end if
        iparenthCount = iparenthCount-1
        if (iparenthCount .lt. 0) then
          call output_line('Mismatched parenthesis in string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        if (sfunctionString(ifunctionIndex-1:ifunctionIndex-1) .eq. '(') then
          call output_line('Empty parentheses in string '//&
              trim(adjustl(sfunctionString))//' !',&
              OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
          call sys_halt()
        end if
        ifunctionIndex = ifunctionIndex+1
        if (ifunctionIndex .gt. ifunctionLength) exit
        c = sfunctionString(ifunctionIndex:ifunctionIndex)
      end do

      ! Now, we have a legal operand: A legal operator or end of string must follow
999   if (ifunctionIndex .gt. ifunctionLength) exit

      ! Check operators
      iopSize = 0
      if (.not.stack_empty(rstack)) then
        call stack_top(rstack, idummy)
        if (c .eq. ',' .and. idummy .eq. iparenthCount) then
          iopSize = 1
        else
          iopSize = isOperator(sfunctionString(ifunctionIndex:))
        end if
      else
        iopSize = isOperator(sfunctionString(ifunctionIndex:))
      end if
      if (iopSize .eq. 0) then
        call output_line('Operator expected in string '//&
              trim(adjustl(sfunctionString))//' !',&
            OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
        call sys_halt()
      end if

      ! Now, we have an operand and an operator: the next loop will check for another
      ! operand (must appear)
      ifunctionIndex = ifunctionIndex+iopSize
    end do
    if (iparenthCount .gt. 0) then
      call output_line('Missing ) in string '//&
              trim(adjustl(sfunctionString))//' !',&
          OU_CLASS_ERROR, OU_MODE_STD,'CheckSyntax')
      call sys_halt()
    end if

    call stack_release(rstack)

  end subroutine CheckSyntax

  ! *****************************************************************************

!<function>

  function isOperator (sfunctionString) result (n)

!<description>
    ! Return 0 if given string is not an operator, else the size of the
    ! operator
!</description>

!<input>
    ! Operator string
    character(LEN=*), intent(in) :: sfunctionString
!</input>

!<result>
    ! Length of operator, 0 if string is no operator
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j,m

    n = 0
    do j = cAdd, cOr
      m = len_trim(Ops(j))
      if (sfunctionString(1:m) .eq. trim(Ops(j))) then
        n = m
        exit
      end if
    end do

  end function isOperator

  ! *****************************************************************************

!<function>

  function MathFunctionIndex (sfunctionString) result (n)

!<description>
    ! Return index of math function beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    character (LEN=*), intent(in) :: sfunctionString
!</input>

!<result>
    ! Index of math function
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j
    integer k
    character (LEN=len(Funcs)) :: sfunctionName

    ! Check all math functions
    n = 0
    do j = cIf, cSign
      k = min(len_trim(Funcs(j)), len(sfunctionString))

      ! Compare lower case letters
      call sys_tolower(sfunctionString(1:k), sfunctionName)
      if (sfunctionName .eq. Funcs(j)) then
        ! Found a matching function
        n = j
        return
      end if
    end do

  end function MathFunctionIndex

  ! *****************************************************************************

!<function>

  function MathFunctionParameters (ifunctionIndex) result (nparameters)

!<description>
    ! Return number of required parameters
!</description>

!<input>
    ! Index of function
    integer(is) :: ifunctionIndex
!</input>

!<result>
    ! Number if required parameters
    integer :: nparameters
!</result>
!</function>

    select case(ifunctionIndex)
    case(cIf)
      nparameters = 3

    case(cMin:cAtan2)
      nparameters = 2

    case(cAbs:cSign)
      nparameters = 1

    case default
      nparameters = 0
      call output_line('Not a function',&
          OU_CLASS_WARNING, OU_MODE_STD,'MathFunctionParameters')
    end select

  end function MathFunctionParameters

  ! *****************************************************************************

!<function>

  function ConstantIndex (sfunctionString) result (n)

!<description>
    ! Return index of predefined constants beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    character (LEN=*), intent(in) :: sfunctionString
!</input>

!<result>
    ! Index of math function
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j
    integer k
    character (LEN=len(CconstantName)) :: sconstantName

    ! Check all predefined constants
    n = 0
    do j = 1, nconstants
      k = min(len_trim(CconstantName(j)), len(sfunctionString(2:)))

      ! Compare lower case letters
      call sys_tolower(sfunctionString(2:k+1), sconstantName)
      if (sconstantName .eq. CconstantName(j)) then
        ! Found a matching constant
        n = j
        return
      end if
    end do

  end function ConstantIndex

  ! *****************************************************************************

!<function>

  function ExpressionIndex (sfunctionString) result (n)

!<description>
    ! Return index of predefined expression beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    character (LEN=*), intent(in) :: sfunctionString
!</input>

!<result>
    ! Index of math function
    integer(is) :: n
!</result>
!</function>

    ! local variables
    integer(is) :: j
    integer k
    character (LEN=len(CexpressionName)) :: sexpressionName

    ! Check all predefined expressions
    n = 0
    do j = 1, nexpressions
      k = min(len_trim(CexpressionName(j)), len(sfunctionString(2:)))

      ! Compare lower case letters
      call sys_tolower(sfunctionString(2:k+1), sexpressionName)

      if (sexpressionName .eq. CexpressionName(j)) then
        ! Found a matching expression
        n = j
        return
      end if
    end do

  end function ExpressionIndex

  ! *****************************************************************************

!<function>

  function VariableIndex (sstring, Svariables, ibegin, inext) result (n)

!<description>
    ! Return index of variable at begin of string sfunctionString (returns 0 if no variable found)
!</description>

!<input>
    ! String
    character (LEN=*), intent(in) :: sstring

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Svariables
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
    integer :: j,ib,in,istringlen

    n = 0
    istringlen = len_trim(sstring)
    if (istringlen .gt. 0) then
      ! Search for first character in str
      do ib = 1, istringlen
        ! When lstr>0 at least 1 char in str
        if (sstring(ib:ib) .ne. ' ') exit
      end do

      ! Search for name terminators
      do in = ib, istringlen
        if (scan(sstring(in:in),'*+-/%^),&|<>=! ') > 0) exit
      end do
      do j = 1, size(Svariables)
        if (sstring(ib:in-1) .eq. Svariables(j)) then
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

  subroutine RemoveSpaces (sfunctionString)

!<description>
    ! Remove Spaces from string, remember positions of characters in
    ! old string
!</description>

!<inputoutput>
    ! String from which spaces should be removed
    character (LEN=*), intent(inout) :: sfunctionString
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: k,istringlen

    istringlen = len_trim(sfunctionString)
    k = 1
    do while (sfunctionString(k:istringlen) .ne. ' ')
      if (sfunctionString(k:k) .eq. ' ') then
        sfunctionString(k:istringlen)  = sfunctionString(k+1:istringlen)//' ' ! Move 1 character to left
        k = k-1
      end if
      k = k+1
    end do

  end subroutine RemoveSpaces

  ! *****************************************************************************

!<subroutine>

  subroutine Replace (ca, cb, sfunctionString)

!<description>
    ! Replace ALL appearances of character set ca in string sfunctionString by character set cb
!</description>

!<input>
    ! Source characters
    character (LEN=*), intent(in) :: ca

    ! Destination characters
    character (LEN=len(ca)), intent(in) :: cb
!</input>

!<inputoutput>
    ! String
    character (LEN=*), intent(inout) :: sfunctionString
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: j,lca

    lca = len(ca)
    do j = 1, len_trim(sfunctionString)-lca+1
      if (sfunctionString(j:j+lca-1) .eq. ca) sfunctionString(j:j+lca-1) = cb
    end do

  end subroutine Replace

  ! *****************************************************************************

!<subroutine>

  subroutine Compile (rcomp, sfunctionString, Svariables)

!<description>
    ! Compile i-th function string into bytecode
!</description>

!<input>
    ! Function string
    character (LEN=*), intent(in) :: sfunctionString

    ! Array with variable names
    character (LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser component
    type (t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>
!</subroutine>

    ! local variables
    integer(is), dimension(:), pointer :: IbyteCode
    real(DP), dimension(:), pointer :: Dimmed
    integer :: ind,isize

    ! (Re-)initialise the bytecode structure (if required)
    if (associated(rcomp%IbyteCode)) deallocate (rcomp%IbyteCode)
    if (associated(rcomp%Dimmed))    deallocate (rcomp%Dimmed)
    rcomp%ibytecodeSize = 0
    rcomp%iimmedSize    = 0
    rcomp%istackSize    = 0
    rcomp%istackPtr     = 0
    rcomp%bisVectorizable = .true.

    ! Neither the stack for the bytecode nor the stack for the
    ! immediate expressions can exceed the size of the function
    ! string. Hence, allocate some initial memory
    isize = FunctionSize(sfunctionString)
    allocate(rcomp%IbyteCode(isize), rcomp%Dimmed(isize))

    ! Compile function string into bytecode
    ind = CompileExpression(rcomp, sfunctionString, 1, Svariables)

    ! Adjust memory size of bytecode stack
    if (rcomp%ibytecodeSize .eq. 0) then
      deallocate(rcomp%IbyteCode)
    else
      allocate(IbyteCode(rcomp%ibytecodeSize))
      IbyteCode = rcomp%IbyteCode(1:rcomp%ibytecodeSize)
      deallocate(rcomp%IbyteCode)
      allocate(rcomp%IbyteCode(rcomp%ibytecodeSize))
      rcomp%IbyteCode = IbyteCode
      deallocate(IbyteCode)
    end if

    ! Adjust memory size of immediate stack
    if (rcomp%iimmedSize .eq. 0) then
      deallocate(rcomp%Dimmed)
    else
      allocate(Dimmed(rcomp%iimmedSize))
      Dimmed = rcomp%Dimmed(1:rcomp%iimmedSize)
      deallocate(rcomp%Dimmed)
      allocate(rcomp%Dimmed(rcomp%iimmedSize))
      rcomp%Dimmed = Dimmed
      deallocate(Dimmed)
    end if

  end subroutine Compile

  ! *****************************************************************************

!<subroutine>

  subroutine incStackPtr (rcomp)

!<description>
    ! Increase stack pointer
!</description>

!<inputoutput>
    type (t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>
!</subroutine>

    rcomp%istackPtr = rcomp%istackPtr+1
    if (rcomp%istackPtr .gt. rcomp%istackSize) rcomp%istackSize = rcomp%istackSize+1

  end subroutine incStackPtr

  ! *****************************************************************************

!<subroutine>

  subroutine AddCompiledByte (rcomp, ibyte)

!<description>
    ! Add compiled byte to bytecode
!</description>

!<input>
    ! Value of byte to be added
    integer(is), intent(in) :: ibyte
!</input>

!<inputoutput>
    type (t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP) :: daux
    integer, dimension(:), allocatable :: p_Irandom1,p_Irandom2
    integer :: iaux

    rcomp%ibytecodeSize = rcomp%ibytecodeSize + 1
    rcomp%IbyteCode(rcomp%ibytecodeSize) = ibyte

    ! Try to optimise the compiled bytecode. Check the bytecode instruction and
    ! compute some values on-the-fly of this is possible
    select case(ibyte)
      !------------------------------------------------------------
      ! Functions
      !------------------------------------------------------------
    case (cAbs)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = abs(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cAcos)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .lt. -1.0_DP .or. &
            rcomp%Dimmed(rcomp%iimmedSize) .gt.  1.0_DP) then
          call output_line('Invalid argument for ACOS!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = acos(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cAsin)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .lt. -1.0_DP .or. &
            rcomp%Dimmed(rcomp%iimmedSize) .gt.  1.0_DP) then
          call output_line('Invalid argument for ASIN!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = asin(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cAtan)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = atan(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cAtan2)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = atan2(rcomp%Dimmed(rcomp%iimmedSize),&
                                                 rcomp%Dimmed(rcomp%iimmedSize-1))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cAcosh)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        daux=rcomp%Dimmed(rcomp%iimmedSize)+sqrt(rcomp%Dimmed(rcomp%iimmedSize)**2-1)
        if (daux .le. 0) then
          call output_line('Invalid argument for ACOSH!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = log(daux)
        call RemoveCompiledByte(rcomp)
      end if

    case (cAnint)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = anint(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cAint)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = aint(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cAsinh)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        daux=rcomp%Dimmed(rcomp%iimmedSize)+sqrt(rcomp%Dimmed(rcomp%iimmedSize)**2-1)
        if (daux .le. 0) then
          call output_line('Invalid argument for ASINH!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = log(daux)
        call RemoveCompiledByte(rcomp)
      end if

    case (cAtanh)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .eq. -1.0_DP) then
          call output_line('Invalid argument for ATANH!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        daux=(1+rcomp%Dimmed(rcomp%iimmedSize))/(1-rcomp%Dimmed(rcomp%iimmedSize))
        if (daux .le. 0._DP) then
          call output_line('Invalid argument for ATANH!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = log(daux)/2.0_DP
        call RemoveCompiledByte(rcomp)
      end if

    case (cCeil)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = ceiling(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cCos)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = cos(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cCosh)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = cosh(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cCot)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        daux=tan(rcomp%Dimmed(rcomp%iimmedSize))
        if (daux .eq. 0.0_DP) then
          call output_line('Invalid argument for COT!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = 1/daux
        call RemoveCompiledByte(rcomp)
      end if

    case (cCsc)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        daux=sin(rcomp%Dimmed(rcomp%iimmedSize))
        if (daux .eq. 0._DP) then
          call output_line('Invalid argument for CSC!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = 1/daux
        call RemoveCompiledByte(rcomp)
      end if

    case (cExp)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = exp(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cFloor)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = floor(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cIf)
      ! No optimization possible

    case (cLog)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .le. 0.0_DP) then
          call output_line('Invalid argument for LOG!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = log(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cLog10)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .le. 0.0_DP) then
          call output_line('Invalid argument for LOG!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = log10(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cMax)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = max(rcomp%Dimmed(rcomp%iimmedSize),&
                                               rcomp%Dimmed(rcomp%iimmedSize-1))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cMin)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = min(rcomp%Dimmed(rcomp%iimmedSize),&
                                               rcomp%Dimmed(rcomp%iimmedSize-1))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cRrand)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        call random_seed (size=iaux)
        allocate (p_Irandom1(iaux))
        allocate (p_Irandom2(iaux))
        call random_seed (get=p_Irandom1)

        p_Irandom2(:) = 0
        p_Irandom2(1) = int(rcomp%Dimmed(rcomp%iimmedSize-1))
        call random_seed (put=p_Irandom2)
        daux = 0.0_DP
        do iaux=1,max(1,int(rcomp%Dimmed(rcomp%iimmedSize)))
          call random_number (daux)
        end do
        rcomp%Dimmed(rcomp%iimmedSize-1) = daux

        call random_seed (put=p_Irandom1)
        deallocate(p_Irandom1,p_Irandom2)

        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cSec)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        daux=cos(rcomp%Dimmed(rcomp%iimmedSize))
        if (daux .eq. 0._DP) then
          call output_line('Invalid argument for SEC!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = 1/daux
        call RemoveCompiledByte(rcomp)
      end if

    case (cSign)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = sign(1._DP,rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cSin)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = sin(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cSinh)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = sinh(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cSqrt)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .lt. 0.0_DP) then
          call output_line('Invalid argument for SQRT!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize) = sqrt(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cTan)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = tan(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cTanh)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = tanh(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

      !------------------------------------------------------------
      ! Misc
      !------------------------------------------------------------
    case (cImmed, cJump)
      ! No optimization needed

      !------------------------------------------------------------
      ! Operators
      !------------------------------------------------------------
    case (cNeg)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = -(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cAdd)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = rcomp%Dimmed(rcomp%iimmedSize-1)+&
                                           rcomp%Dimmed(rcomp%iimmedSize)
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cSub)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = rcomp%Dimmed(rcomp%iimmedSize-1)-&
                                           rcomp%Dimmed(rcomp%iimmedSize)
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cMul)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = rcomp%Dimmed(rcomp%iimmedSize-1)*&
                                           rcomp%Dimmed(rcomp%iimmedSize)
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cDiv)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .eq. 0.0_DP) then
          call output_line('Invalid argument for DIV!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize-1) = rcomp%Dimmed(rcomp%iimmedSize-1)/&
                                           rcomp%Dimmed(rcomp%iimmedSize)
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cMod)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        if (rcomp%Dimmed(rcomp%iimmedSize) .eq. 0.0_DP) then
          call output_line('Invalid argument for MOD!',&
              OU_CLASS_ERROR, OU_MODE_STD,'AddCompiledByte')
          call sys_halt()
        end if
        rcomp%Dimmed(rcomp%iimmedSize-1) = mod(rcomp%Dimmed(rcomp%iimmedSize-1),&
                                               rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cPow)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = rcomp%Dimmed(rcomp%iimmedSize-1)**rcomp%Dimmed(rcomp%iimmedSize)
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cEqual)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(rcomp%Dimmed(rcomp%iimmedSize-1) .eq.&
                                                      rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cNEqual)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(rcomp%Dimmed(rcomp%iimmedSize-1) .ne.&
                                                      rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cLess)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(rcomp%Dimmed(rcomp%iimmedSize-1) .lt.&
                                                      rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cLessOrEq)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(rcomp%Dimmed(rcomp%iimmedSize-1) .le.&
                                                      rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cGreater)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(rcomp%Dimmed(rcomp%iimmedSize-1) .gt.&
                                                      rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cGreaterOrEq)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(rcomp%Dimmed(rcomp%iimmedSize-1) .ge.&
                                                      rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cAnd)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(DbleToLogc(rcomp%Dimmed(rcomp%iimmedSize-1)) .and.&
                                                      DbleToLogc(rcomp%Dimmed(rcomp%iimmedSize)))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cOr)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed .and.&
          rcomp%IbyteCode(rcomp%ibytecodeSize-2) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize-1) = LogcToDble(DbleToLogc(rcomp%Dimmed(rcomp%iimmedSize-1)) .or.&
                                                      DbleToLogc(rcomp%Dimmed(rcomp%iimmedSize)))
        call RemoveCompiledImmediate(rcomp)
        call RemoveCompiledByte(rcomp)
        call RemoveCompiledByte(rcomp)
      end if

    case (cNot)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = LogcToDble(.not.DbleToLogc(rcomp%Dimmed(rcomp%iimmedSize)))
        call RemoveCompiledByte(rcomp)
      end if

      !------------------------------------------------------------
      ! Degrees-radians conversion
      !------------------------------------------------------------
    case (cDeg)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = RadToDeg(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if

    case (cRad)
      if (rcomp%IbyteCode(rcomp%ibytecodeSize-1) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = DegToRad(rcomp%Dimmed(rcomp%iimmedSize))
        call RemoveCompiledByte(rcomp)
      end if
    end select

  end subroutine AddCompiledByte

  ! *****************************************************************************

!<subroutine>

  subroutine RemoveCompiledByte (rcomp)

!<description>
    ! Remove last compiled byte from bytecode
!</description>

!<inputoutput>
    type (t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>
!</subroutine>

    rcomp%IbyteCode(rcomp%ibytecodeSize) = 0
    rcomp%ibytecodeSize = rcomp%ibytecodeSize - 1

  end subroutine RemoveCompiledByte

  ! *****************************************************************************

!<subroutine>

  subroutine AddImmediate (rcomp, immediate)

!<description>
    ! Add immediate
!</description>

!<input>
    ! Value of byte to be added
    real(DP), intent(in) :: immediate
!</input>

!<inputoutput>
    type (t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>
!</subroutine>

    rcomp%iimmedSize = rcomp%iimmedSize + 1
    rcomp%Dimmed(rcomp%iimmedSize) = immediate

  end subroutine AddImmediate

  ! *****************************************************************************

!<subroutine>

  subroutine RemoveCompiledImmediate (rcomp)

!<description>
    ! Remove last compiled immediate from immediate stack
!</description>

!<inputoutput>
    type (t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>
!</subroutine>

    rcomp%Dimmed(rcomp%iimmedSize) = 0
    rcomp%iimmedSize = rcomp%iimmedSize - 1

  end subroutine RemoveCompiledImmediate

  ! *****************************************************************************

!<subroutine>

  subroutine AddFunctionOpcode (rcomp, iopcode)

!<description>
    ! Add compiled byte to bytecode
!</description>

!<input>
    ! Value of opcode to be added
    integer(is), intent(in) :: iopcode
!</input>

!<inputoutput>
    type (t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>
!</subroutine>

    if (rcomp%buseDegreeConversion) then
      select case(iopcode)
      case(cCos, Ccosh, cCot, cCsc, cSec, cSin, cSinh, cTan, cTanh)
        call AddCompiledByte(rcomp, cRad)
      end select
    end if

    call AddCompiledByte(rcomp, iopcode)

    if (rcomp%buseDegreeConversion) then
      select case(iopcode)
      case(cAcos, cAcosh, cAsinh, cAtanh, cAsin, cAtan, cAtan2)
        call AddCompiledByte(rcomp, cDeg)
      end select
    end if

  end subroutine AddFunctionOpcode

  ! *****************************************************************************

!<function>

  function RealNum (sfunctionString, ibegin, inext, berror) result (dresult)

!<description>
    ! Get real number from string
    ! Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
!</description>

!<input>
    ! String
    character (LEN=*), intent(in) :: sfunctionString
!</input>

!<output>
    ! OPTIONAL: Start position of real number
    integer, optional, intent(out) :: ibegin

    ! OPTIONAL: 1st character after real number
    integer, optional, intent(out) :: inext

    ! OPTIONAL: Error flag
    logical, optional, intent(out) :: berror
!</output>

!<result>
    ! Real number
    real(DP) :: dresult
!</result>
!</function>

    ! local variables
    integer :: ib,in,istat
    logical :: Bflag,               & ! .T. at begin of number in string
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
    do while (in .le. len_trim(sfunctionString))
      select case (sfunctionString(in:in))
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
      case default
        exit ! CALL sys_halt() at all other characters
      end select
      in = in+1
    end do
    err = (ib .gt. in-1) .or. (.not.DInMan) .or.&
          ((Eflag.or.InExp).and..not.DInExp)
    if (err) then
      dresult = 0.0_DP
    else
      read(sfunctionString(ib:in-1),*, IOSTAT=istat) dresult
      err = istat .ne. 0
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
    if (present(berror)) berror = err

  end function RealNum

  ! *****************************************************************************

!<function>

  recursive function FunctionSize (sfunctionString) result (isize)

!<description>
    ! Return the size of the total function including external expressions
!</description>

!<input>
    ! Function string
    character (LEN=*), intent(in) :: sfunctionString
!</input>

!<result>
    ! Size of function string
    integer :: isize
!</result>
!</function>

    ! local variables
    integer :: ind,n
    character(LEN=1) :: c

    ! Determine size of given expression
    isize = len_trim(sfunctionString)

    ! "Parse" string for externally defined expressions
    do ind = 1, isize
      c = sfunctionString(ind:ind)
      if (c .eq. '@') then
        n = ExpressionIndex (sfunctionString(ind:))
        if (n .eq. 0) then
          call output_line('Invalid expression!',&
              OU_CLASS_ERROR, OU_MODE_STD,'FunctionSize')
          call sys_halt()
        end if
        isize = isize+FunctionSize(CexpressionString(n))
      end if
    end do

  end function FunctionSize

  ! *****************************************************************************

!<function>

  recursive function CompileExpression(rcomp, sfunctionString, ind, Svariables,&
                                       bstopAtComma) result(ind2)

!<description>
    ! Compiles ','
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables

    ! OPTIONAL: stop at comma
    logical, intent(in), optional :: bstopAtComma
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileOr(rcomp, sfunctionString, ind, Svariables)
    if(ind2 > ifunctionLength) return

    if (present(bstopAtComma)) then
      if (bstopAtComma) return
    end if

    do while (sfunctionString(ind2:ind2) .eq. ',')
      ind2 = CompileOr(rcomp, sfunctionString, ind2+1, Svariables)
      if (ind2 .gt. ifunctionLength) return
    end do

  end function CompileExpression

  ! *****************************************************************************

!<function>

  recursive function CompileOr(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles '|'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileAnd(rcomp, sfunctionString, ind, Svariables)
    if (ind2 .gt. ifunctionLength) return

    do while(sfunctionString(ind2:ind2) .eq. '|')
      ind2 = CompileAnd(rcomp, sfunctionString, ind2+1, Svariables)

      call AddCompiledByte(rcomp, cOr)
      rcomp%istackPtr = rcomp%istackPtr-1
      if (ind2 .gt. ifunctionLength) return
    end do

  end function CompileOr

  ! *****************************************************************************

!<function>

  recursive function CompileAnd(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles '&'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileComparison(rcomp, sfunctionString, ind, Svariables)
    if (ind2 .gt. ifunctionLength) return

    do while(sfunctionString(ind2:ind2) .eq. '&')
      ind2 = CompileComparison(rcomp, sfunctionString, ind2+1, Svariables)

      call AddCompiledByte(rcomp, cAnd)
      rcomp%istackPtr = rcomp%istackPtr-1
      if (ind2 .gt. ifunctionLength) return
    end do

  end function CompileAnd

  ! *****************************************************************************

!<function>

  recursive function CompileComparison(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles '=', '<' and '>'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer(is) :: iopSize
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileAddition(rcomp, sfunctionString, ind, Svariables)
    if (ind2 .gt. ifunctionLength) return

    c=sfunctionString(ind2:ind2)
    do while(c .eq. '=' .or. c .eq. '<' .or. c .eq. '>' .or. c .eq. '!')
      iopSize = merge(2, 1, sfunctionString(ind2+1:ind2+1) .eq. '=')
      ind2 = CompileAddition(rcomp, sfunctionString, ind2+iopSize, Svariables)

      select case(c)
      case('=')
        call AddCompiledByte(rcomp, cEqual)

      case('<')
        call AddCompiledByte(rcomp, merge(cLess, cLessOrEq, iopSize .eq. 1))

      case('>')
        call AddCompiledByte(rcomp, merge(cGreater, cGreaterOrEq, iopSize .eq. 1))

      case('!')
        call AddCompiledByte(rcomp, cNEqual)
      end select
      rcomp%istackPtr = rcomp%istackPtr-1

      if (ind2 .gt. ifunctionLength) return
      c=sfunctionString(ind2:ind2)
    end do

  end function CompileComparison

  ! *****************************************************************************

!<function>

  recursive function CompileAddition(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles '+' and '-'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileMult(rcomp, sfunctionString, ind, Svariables)
    if (ind2 .gt. ifunctionLength) return

    c=sfunctionString(ind2:ind2)
    do while(c .eq. '+' .or. c .eq. '-')
      ind2 = CompileMult(rcomp, sfunctionString, ind2+1, Svariables)

      call AddCompiledByte(rcomp, merge(cAdd, cSub, c .eq. '+'))
      rcomp%istackPtr = rcomp%istackPtr-1

      if (ind2 .gt. ifunctionLength) return
      c=sfunctionString(ind2:ind2)
    end do

  end function CompileAddition

  ! *****************************************************************************

!<function>

  recursive function CompileMult(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles '*', '/' and '%'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileUnaryMinus(rcomp, sfunctionString, ind, Svariables)
    if (ind2 .gt. ifunctionLength) return

    c=sfunctionString(ind2:ind2)
    do while(c .eq. '*' .or. c .eq. '/' .or. c .eq. '%')
      ind2 = CompileUnaryMinus(rcomp, sfunctionString, ind2+1, Svariables)

      select case(c)
      case('*')
        call AddCompiledByte(rcomp, cMul)

      case('/')
        call AddCompiledByte(rcomp, cDiv)

      case('%')
        call AddCompiledByte(rcomp, cMod)

      end select
      rcomp%istackPtr = rcomp%istackPtr-1

      if (ind2 .gt. ifunctionLength) return
      c=sfunctionString(ind2:ind2)
    end do

  end function CompileMult

  ! *****************************************************************************

!<function>

  recursive function CompileUnaryMinus(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles unary '-'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    c=sfunctionString(ind:ind)
    if (c .eq. '-' .or. c .eq. '!') then
      ind2 = ind+1
      if (ind2 .gt. ifunctionLength) return
      ind2 = CompilePow(rcomp, sfunctionString, ind2, Svariables)

      ! If we are negating a constant, negate the constant itself
      if (c .eq. '-' .and. rcomp%IbyteCode(rcomp%ibytecodeSize) .eq. cImmed) then
        rcomp%Dimmed(rcomp%iimmedSize) = -rcomp%Dimmed(rcomp%iimmedSize)

        ! If we are negating a negation, we can remove both
      elseif (c .eq. '-' .and. rcomp%IbyteCode(rcomp%ibytecodeSize) .eq. cNeg) then
        call RemoveCompiledByte(rcomp)

      else
        call AddCompiledByte(rcomp, merge(cNeg, cNot, c .eq. '-'))

      end if
      return
    end if

    ind2 = CompilePow(rcomp, sfunctionString, ind, Svariables)

  end function CompileUnaryMinus

  ! *****************************************************************************

!<function>

  recursive function CompilePow(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles '^'
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: ifunctionLength

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileElement(rcomp, sfunctionString, ind, Svariables)
    if (ind2 .gt. ifunctionLength) return

    do while(sfunctionString(ind2:ind2) .eq. '^')
      ind2 = CompileUnaryMinus(rcomp, sfunctionString, ind2+1, Svariables)

      call AddCompiledByte(rcomp, cPow)
      rcomp%istackPtr = rcomp%istackPtr-1
      if (ind2 .gt. ifunctionLength) return
    end do

  end function CompilePow

  ! *****************************************************************************

!<function>

  recursive function CompileElement(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles element
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    character(LEN=1) :: c
    real(DP) :: dnumber
    integer(is) :: n
    integer :: ind1,ind0,ib,in,nparams
    logical :: berror

    ind1=ind; c=sfunctionString(ind1:ind1)
    if (c .eq. '(') then
      ind1 = CompileExpression(rcomp, sfunctionString, ind1+1, Svariables, .false.)
      ind2 = ind1+1   ! sfunctionString(ind1:ind1) is ')'
      return
    end if

    ! Check for numbers
    if (scan(c,'0123456789,') > 0) then
      dnumber = RealNum (sfunctionString(ind1:), ib, in, berror)
      if (berror) then
        call output_line('Invalid number format!',&
            OU_CLASS_ERROR, OU_MODE_STD,'CompileElement')
        call sys_halt()
      end if
      call AddImmediate(rcomp, dnumber)
      call AddCompiledByte(rcomp, cImmed)
      call incstackPtr(rcomp)
      ind2 = ind1+in-1
      return

    else
      ! Then string must be function, variable or constant

      ! Check for mathematical functions
      n = MathFunctionIndex(sfunctionString(ind1:))
      if (n .gt. 0) then
        ind2 = ind1+len_trim(Funcs(n))

        ! Check for IF-THEN-ELSE
        if (n .eq. cIf) then
          ind2 = CompileIf(rcomp, sfunctionString, ind2+1, Svariables)
          ! IF-THEN-ELSE cannot be vectorised, note that!
          rcomp%bisVectorizable = .false.
          return
        end if

        nparams = MathFunctionParameters(n)
        ind2 = CompileFunctionParameters(rcomp, sfunctionString, ind2+1, Svariables, nparams)
        call AddFunctionOpcode(rcomp, n)
        return
      end if

      ! Check for predefined constant
      n = ConstantIndex(sfunctionString(ind1:))
      if (n .gt. 0) then
        ind2 = ind1+len_trim(CconstantName(n))+1
        call AddImmediate(rcomp, DconstantValue(n))
        call AddCompiledByte(rcomp, cImmed)
        call incStackPtr(rcomp)
        return
      end if

      ! Check for predefined expressions
      n = ExpressionIndex(sfunctionString(ind1:))
      if (n .gt. 0) then
        ind2 = ind1+len_trim(CexpressionName(n))+1

        ! Recursively compile the given expression
        ind0 = CompileExpression(rcomp, CexpressionString(n), 1, Svariables)

        ! Proceed with compilation of mathematical function Func afterwards
        return
      end if

      ! Check for variables
      n = VariableIndex(sfunctionString(ind1:), Svariables, ib, in)
      if (n > 0) n = VarBegin+n-1
      call AddCompiledByte(rcomp, n)
      call incStackPtr(rcomp)
      ind2 = ind1+in-1
      return
    end if

    call output_line('An unexpected error occured!',&
        OU_CLASS_ERROR, OU_MODE_STD,'CompileElement')
    call sys_halt()

  end function CompileElement

  ! *****************************************************************************

!<function>

  recursive function CompileFunctionParameters(rcomp, sfunctionString, ind, Svariables,&
                                               nparams) result(ind2)

!<description>
    ! Compiles function parameters
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables

    ! Number of required parameters
    integer, intent(in) :: nparams
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: iStackPtr

    ind2 = ind
    if (nparams .gt. 0) then

      iStackPtr = rcomp%istackPtr
      ind2 = CompileExpression(rcomp, sfunctionString, ind, Svariables, .false.)

      if (rcomp%istackPtr .ne. iStackPtr+nparams) then
        call output_line('Illegal number of parameters to function!',&
            OU_CLASS_ERROR, OU_MODE_STD,'CompileFunctionParameters')
        call sys_halt()
      end if

      rcomp%istackPtr = rcomp%istackPtr-(nparams-1)

    else

      call incStackPtr(rcomp)

    end if
    ind2=ind2+1

  end function CompileFunctionParameters

  ! *****************************************************************************

!<function>

  recursive function CompileIf(rcomp, sfunctionString, ind, Svariables) result(ind2)

!<description>
    ! Compiles if()
!</description>

!<input>
    ! Function substring
    character(LEN=*), intent(in) :: sfunctionString

    ! Begin position substring
    integer, intent(in) :: ind

    ! Array with variable names
    character(LEN=*), dimension(:), intent(in) :: Svariables
!</input>

!<inputoutput>
    ! Function parser
    type(t_fparserComponent), intent(inout) :: rcomp
!</inputoutput>

!<result>
    integer :: ind2
!</result>
!</function>

    ! local variables
    integer :: ifunctionLength,curibytecodeSize,curibytecodeSize2,curiimmedSize2

    ifunctionLength = len_trim(sfunctionString)

    ind2 = CompileExpression(rcomp, sfunctionString, ind, Svariables, .true.) ! Condition branch
    if (ind2 .gt. ifunctionLength) return

    if (sfunctionString(ind2:ind2) .ne. ',') then
      call output_line('Illegal number of parameters to function!',&
          OU_CLASS_ERROR, OU_MODE_STD,'CompileIf')
      call sys_halt()
    end if
    call AddCompiledByte(rcomp, cIf)
    curibytecodeSize = rcomp%ibytecodeSize
    call AddCompiledByte(rcomp, 0_is) ! Jump index will be set below
    call AddCompiledByte(rcomp, 0_is) ! Immed jump index will be set below
    rcomp%istackPtr = rcomp%istackPtr-1

    ind2 = CompileExpression(rcomp, sfunctionString, ind2+1, Svariables, .true.) ! Then branch
    if (ind2 .gt. ifunctionLength) return

    if (sfunctionString(ind2:ind2) .ne. ',') then
      call output_line('Illegal number of parameters to function!',&
          OU_CLASS_ERROR, OU_MODE_STD,'CompileIf')
      call sys_halt()
    end if
    call AddCompiledByte(rcomp, cJump)
    curibytecodeSize2 = rcomp%ibytecodeSize
    curiimmedSize2 = rcomp%iimmedSize
    call AddCompiledByte(rcomp, 0_is) ! Jump index will be set below
    call AddCompiledByte(rcomp, 0_is) ! Immed jump index will be set below
    rcomp%istackPtr = rcomp%istackPtr-1

    ind2 = CompileExpression(rcomp, sfunctionString, ind2+1, Svariables, .true.) ! Else branch
    if (ind2 .gt. ifunctionLength) return

    if (sfunctionString(ind2:ind2) .ne. ')') then
      call output_line('Illegal number of parameters to function!',&
          OU_CLASS_ERROR, OU_MODE_STD,'CompileIf')
      call sys_halt()
    end if

    ! Set jump indices
    if (associated(rcomp%IbyteCode)) then
      rcomp%IbyteCode(curibytecodeSize+1)  = curibytecodeSize2+2
      rcomp%IbyteCode(curibytecodeSize+2)  = curiimmedSize2+1
      rcomp%IbyteCode(curibytecodeSize2+1) = rcomp%ibytecodeSize
      rcomp%IbyteCode(curibytecodeSize2+2) = rcomp%iimmedSize+1
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
    real(DP), intent(in) :: d
!</input>

!<result>
    ! Logical variable
    logical :: l
!</result>
!</function>

    l = (abs(1-d) .le. 1e-12)

  end function DbleTOLogc

  ! *****************************************************************************

!<function>

  elemental function LogcToDble(l) result(d)

!<description>
    ! This function transforms a Logical into a Double
!</description>

!<input>
    ! Logical variable
    logical, intent(in) :: l
!</input>

!<result>
    ! Double variable
    real(DP) :: d
!</result>
!</function>

    d = merge(1._DP, 0._DP, l)

  end function LogcToDble

  ! *****************************************************************************

!<function>

  elemental function DegToRad(d) result(r)

!<description>
    ! This function converts DEG to RAD
!</description>

!<input>
    ! DEG
    real(DP), intent(in) :: d
!</input>

!<result>
    ! RAD
    real(DP) :: r
!</result>
!</function>

    r = d * (3.141592653589793115997963468544185161590576171875_DP / 180._DP)
  end function DegToRad

  ! *****************************************************************************

!<function>

  elemental function RadToDeg(r) result(d)

!<description>
    ! This function converts RAD to DEG
!</description>

!<input>
    ! RAD
    real(DP), intent(in) :: r
!</input>

!<result>
    ! DEG
    real(DP) :: d
!</result>
!</function>

    d = r * (180._DP / 3.141592653589793115997963468544185161590576171875_DP)

  end function RadToDeg

  ! *****************************************************************************

!<subroutine>

  subroutine evalFunctionScalar (rcomp, Dstack, Dvalue, EvalErrType, dresult)

!<description>
    ! Evaluate bytecode for the values passed in array Val(:).
!</description>

!<input>
    ! Component of function parser
    type(t_fparserComponent), intent(in) :: rcomp

    ! Variable values
    real(DP), dimension(:), intent(in) :: Dvalue
!</input>

!<inputoutput>
    ! Stack memory
    real(DP), dimension(:), intent(inout) :: Dstack
!</inputoutput>

!<output>
    ! Error code for function evaluation
    integer, intent(out) :: EvalErrType

    ! Evaluated function
    real(DP), intent(out) :: dresult
!</output>
!</subroutine>

    ! local variables
    integer  :: iinstPtr,istackPtr,idataPtr
    integer  :: ijumpAddr,iimmedAddr
    real(DP) :: daux
    integer :: iaux
    integer, dimension(:), allocatable :: p_Irandom1,p_Irandom2

    ! Initialization
    idataPtr  = 1
    istackPtr = 0
    iinstPtr  = 0

    ! Repeat until complete bytecode has been processed
    do while(iinstPtr .lt. rcomp%ibytecodeSize)
      iinstPtr = iinstPtr+1

      ! What kind of bytecode are we?
      select case (rcomp%IbyteCode(iinstPtr))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      case (cAbs)
        Dstack(istackPtr) = abs(Dstack(istackPtr))

      case (cAcos)
        if ((Dstack(istackPtr) < -1._DP) .or. &
            (Dstack(istackPtr) >  1._DP)) then
          EvalErrType = 4
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr) = acos(Dstack(istackPtr))

      case (cAsin)
        if ((Dstack(istackPtr) < -1._DP) .or. &
            (Dstack(istackPtr) >  1._DP)) then
          EvalErrType = 4
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr) = asin(Dstack(istackPtr))

      case (cAtan)
        Dstack(istackPtr) = atan(Dstack(istackPtr))

      case (cAtan2)
        Dstack(istackPtr-1) = atan2(Dstack(istackPtr -1), Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cAcosh)
        daux = Dstack(istackPtr)+sqrt(Dstack(istackPtr)**2-1)
        if (daux .le. 0._DP) then
          EvalErrType = 5
          dresult     = 0._DP
          return
        end if
        Dstack(istackPtr) = log(daux)

      case (cAnint)
        Dstack(istackPtr) = anint(Dstack(istackPtr))

      case (cAint)
        Dstack(istackPtr) = aint(Dstack(istackPtr))

      case (cAsinh)
        daux = Dstack(istackPtr)+sqrt(Dstack(istackPtr)**2+1)
        if (daux .le. 0._DP) then
          EvalErrType = 5
          dresult     = 0._DP
          return
        end if
        Dstack(istackPtr) = log(daux)

      case (cAtanh)
        if (Dstack(istackPtr) .eq. -1._DP) then
          EvalErrType = 6
          dresult     = 0._DP
          return
        end if
        daux = (1+Dstack(istackPtr))/(1-Dstack(istackPtr))
        if (daux .le. 0._DP) then
          EvalErrType = 3
          dresult     = 0._DP
          return
        end if
        Dstack(istackPtr) = log(daux)/2._DP

      case (cCeil)
        Dstack(istackPtr) = ceiling(Dstack(istackPtr))

      case (cCos)
        Dstack(istackPtr) = cos(Dstack(istackPtr))

      case (cCosh)
        Dstack(istackPtr) = cosh(Dstack(istackPtr))

      case (cCot)
        daux = tan(Dstack(istackPtr))
        if (daux .eq. 0._DP) then
          EvalErrType = 1
          dresult     = 0._DP
          return
        end if
        Dstack(istackPtr) = 1._DP/daux

      case (cCsc)
        daux = sin(Dstack(istackPtr))
        if (daux .eq. 0._DP) then
          EvalErrType = 1
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr) = 1._DP/daux

      case (cExp)
        Dstack(istackPtr) = exp(Dstack(istackPtr))

      case (cFloor)
        Dstack(istackPtr) = floor(Dstack(istackPtr))

      case (cIf)
        iinstPtr = iinstPtr+1;   ijumpAddr  = rcomp%IbyteCode(iinstPtr)
        iinstPtr = iinstPtr+1;   iimmedAddr = rcomp%IbyteCode(iinstPtr)
        if (.not.DbleToLogc(Dstack(istackPtr))) then
          iinstPtr = ijumpAddr
          idataPtr = iimmedAddr
        end if
        istackPtr = istackPtr-1

      case (cLog)
        if (Dstack(istackPtr) .le. 0._DP) then
          EvalErrType = 3
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr) = log(Dstack(istackPtr))

      case (cLog10)
        if (Dstack(istackPtr) .le. 0._DP) then
          EvalErrType = 3
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr) = log10(Dstack(istackPtr))

      case (cMax)
        Dstack(istackPtr-1) = max(Dstack(istackPtr-1), Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cMin)
        Dstack(istackPtr-1) = min(Dstack(istackPtr-1), Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cRrand)
        call random_seed (size=iaux)
        allocate (p_Irandom1(iaux))
        allocate (p_Irandom2(iaux))
        call random_seed (get=p_Irandom1)

        p_Irandom2(:) = 0
        p_Irandom2(1) = int(Dstack(istackPtr-1))
        call random_seed (put=p_Irandom2)

        daux = 0.0_DP
        do iaux=1,max(1,int(Dstack(istackPtr)))
          call random_number (daux)
        end do
        Dstack(istackPtr-1) = daux

        call random_seed (put=p_Irandom1)
        deallocate(p_Irandom1,p_Irandom2)

        istackPtr = istackPtr-1

      case (cSec)
        daux = cos(Dstack(istackPtr))
        if (daux .eq. 0._DP) then
          EvalErrType = 1
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr) = 1._DP/daux

      case (cSign)
        Dstack(istackPtr) = sign(1._DP,Dstack(istackPtr))

      case (cSin)
        Dstack(istackPtr) = sin(Dstack(istackPtr))

      case(cSinh)
        Dstack(istackPtr) = sinh(Dstack(istackPtr))

      case(cSqrt)
        if (Dstack(istackPtr) .lt. 0._DP) then
          EvalErrType = 3
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr) = sqrt(Dstack(istackPtr))

      case (cTan)
        Dstack(istackPtr) = tan(Dstack(istackPtr))

      case (cTanh)
        Dstack(istackPtr) = tanh(Dstack(istackPtr))

      !------------------------------------------------------------
      ! Misc
      !------------------------------------------------------------
      case (cImmed)
        istackPtr         = istackPtr+1
        Dstack(istackPtr) = rcomp%Dimmed(idataPtr)
        idataPtr          = idataPtr+1

      case (cJump)
        idataPtr = rcomp%IbyteCode(iinstPtr+2)
        iinstPtr = rcomp%IbyteCode(iinstPtr+1)

      !------------------------------------------------------------
      ! Operators
      !------------------------------------------------------------
      case (cNeg)
        Dstack(istackPtr) = -Dstack(istackPtr)

      case (cAdd)
        Dstack(istackPtr-1) = Dstack(istackPtr-1)+Dstack(istackPtr)
        istackPtr = istackPtr-1

      case (cSub)
        Dstack(istackPtr-1) = Dstack(istackPtr-1)-Dstack(istackPtr)
        istackPtr = istackPtr-1

      case (cMul)
        Dstack(istackPtr-1) = Dstack(istackPtr-1)*Dstack(istackPtr)
        istackPtr = istackPtr-1

      case (cDiv)
        if (Dstack(istackPtr) .eq. 0._DP) then
          EvalErrType = 1
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr-1) = Dstack(istackPtr-1)/Dstack(istackPtr)
        istackPtr = istackPtr-1

      case (cMod)
        if (Dstack(istackPtr) .eq. 0._DP) then
          EvalErrType = 1
          dresult     = 0._DP
          return
        endif
        Dstack(istackPtr-1) = mod(Dstack(istackPtr-1), Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cPow)
        Dstack(istackPtr-1) = Dstack(istackPtr-1)**Dstack(istackPtr)
        istackPtr = istackPtr-1

      case (cEqual)
        Dstack(istackPtr-1) = LogcToDble(Dstack(istackPtr-1) .eq. Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cNEqual)
        Dstack(istackPtr-1) = LogcToDble(Dstack(istackPtr-1) .ne. Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cLess)
        Dstack(istackPtr-1) = LogcToDble(Dstack(istackPtr-1) .lt. Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cLessOrEq)
        Dstack(istackPtr-1) = LogcToDble(Dstack(istackPtr-1) .le. Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cGreater)
        Dstack(istackPtr-1) = LogcToDble(Dstack(istackPtr-1) .gt. Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cGreaterOrEq)
        Dstack(istackPtr-1) = LogcToDble(Dstack(istackPtr-1) .ge. Dstack(istackPtr))
        istackPtr = istackPtr-1

      case (cAnd)
        Dstack(istackPtr-1) = LogcToDble(DbleToLogc(Dstack(istackPtr-1)) .and. &
                                         DbleToLogc(Dstack(istackPtr)) )
        istackPtr = istackPtr-1

      case (cOr)
        Dstack(istackPtr-1) = LogcToDble(DbleToLogc(Dstack(istackPtr-1)) .or. &
                                         DbleToLogc(Dstack(istackPtr)) )
        istackPtr = istackPtr-1

      case (cNot)
        Dstack(istackPtr) = LogcToDble( .not. DbleToLogc(Dstack(istackPtr)) )

      !------------------------------------------------------------
      ! Degrees-radians conversion
      !------------------------------------------------------------
      case (cDeg)
        Dstack(istackPtr) = RadToDeg(Dstack(istackPtr))

      case (cRad)
        Dstack(istackPtr) = DegToRad(Dstack(istackPtr))

      case default
        istackPtr = istackPtr+1
        Dstack(istackPtr) = DValue(rcomp%IbyteCode(iinstPtr)-VarBegin+1)
      end select
    end do

    EvalErrType = 0
    dresult = Dstack(istackPtr)

  end subroutine evalFunctionScalar

  ! *****************************************************************************

!<subroutine>

  subroutine evalFunctionBlock (rcomp, iblockSize, Dstack, DvalueBlock, idim,&
                                EvalErrType, Dresult, DvalueScalar)

!<description>
    ! Evaluate bytecode for an array of values passed in DvalueBlock(:,:).
!</description>

!<input>
    ! Component of function parser
    type(t_fparserComponent), intent(in) :: rcomp

    ! Variable values
    real(DP), dimension(:,:), intent(in) :: DvalueBlock

    ! Size of the vector block
    integer, intent(in) :: iblockSize

    ! Orientation of the stored values
    ! idim =1 : DvalueBlock is organised as (x1:xN),(y1:yN),...
    ! idim =2 : DvalueBlock is organised as (x1,y1),(x2,y2),...,(xN,yN)
    integer, intent(in) :: idim

    ! Vector of scalar variable values
    real(DP), dimension(:), intent(in), optional :: DvalueScalar
!</input>

!<inputoutput>
    ! Stack memory
    real(DP), dimension(:,:), intent(inout) :: Dstack
!</inputoutput>

!<output>
    ! Error code for function evaluation
    integer, intent(out) :: EvalErrType

    ! Evaluated function
    real(DP), dimension(:), intent(out) :: Dresult
!</output>
!</subroutine>

    ! local variables
    integer  :: iinstPtr,idataPtr,istackPtr,iblock,istartValueScalar,iVariable
    real(DP) :: daux
    integer :: iaux
    integer, dimension(:), allocatable :: p_Irandom1,p_Irandom2

    ! Initialization
    idataPtr  = 1
    istackPtr = 0
    iinstPtr  = 0

    ! This is tricky. istartValueScalar indicates the number of the first
    ! variable which is passed as scalar. Hence, if the optional parameter
    ! DvalueScalar is missing, then istartValueScalar pointers to SIZE(DvalueBlock)+1.
    ! Obviously, no variable beyond this value is addressed.
    istartValueScalar = size(DvalueBlock,3-idim)+1

    ! Repeat until complete bytecode has been processed
    do while(iinstPtr .lt. rcomp%ibytecodeSize)
      iinstPtr = iinstPtr+1

      select case (rcomp%IbyteCode(iinstPtr))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      case (cAbs)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = abs(Dstack(iblock, istackPtr))
        end do


      case (cAcos)
        do iblock = 1, iblockSize
          if ((Dstack(iblock, istackPtr) .lt. -1._DP) .or.&
              (Dstack(iblock, istackPtr) .gt.  1._DP)) then
            EvalErrType = 4
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = acos(Dstack(iblock, istackPtr))
          end if
        end do


      case (cAsin)
        do iblock = 1, iblockSize
          if ((Dstack(iblock, istackPtr) .lt. -1._DP) .or.&
              (Dstack(iblock, istackPtr) .gt.  1._DP)) then
            EvalErrType = 4
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = asin(Dstack(iblock, istackPtr))
          end if
        end do


      case (cAtan)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = atan(Dstack(iblock, istackPtr))
        end do


      case (cAtan2)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = atan2(Dstack(iblock, istackPtr -1),&
                                               Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cAcosh)
        do iblock = 1, iblockSize
          daux = Dstack(iblock, istackPtr)+sqrt(Dstack(iblock, istackPtr)**2-1)
          if (daux .le. 0._DP) then
            EvalErrType = 5
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = log(daux)
          end if
        end do


      case (cAnint)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = anint(Dstack(iblock, istackPtr))
        end do


      case (cAint)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = aint(Dstack(iblock, istackPtr))
        end do


      case (cAsinh)
        do iblock = 1, iblockSize
          daux = Dstack(iblock, istackPtr)+sqrt(Dstack(iblock, istackPtr)**2+1)
          if (daux .le. 0._DP) then
            EvalErrType = 5
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = log(daux)
          end if
        end do


      case (cAtanh)
        do iblock = 1, iblockSize
          if (Dstack(iblock, istackPtr) .eq. -1._DP) then
            EvalErrType = 6
            Dresult(iblock) = 0._DP
          end if
          daux = (1+Dstack(iblock, istackPtr))/(1-Dstack(iblock, istackPtr))
          if (daux .le. 0._DP) then
            EvalErrType = 3
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = log(daux)/2._DP
          end if
        end do


      case (cCeil)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = ceiling(Dstack(iblock, istackPtr))
        end do


      case (cCos)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = cos(Dstack(iblock, istackPtr))
        end do


      case (cCosh)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = cosh(Dstack(iblock, istackPtr))
        end do


      case (cCot)
        do iblock = 1, iblockSize
          daux=tan(Dstack(iblock, istackPtr))
          if (daux .eq. 0) then
            EvalErrType = 1
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = 1._DP/daux
          end if
        end do


      case (cCsc)
        do iblock = 1, iblockSize
          daux=sin(Dstack(iblock, istackPtr))
          if (daux.eq.0._DP) then
            EvalErrType = 1
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = 1._DP/daux
          end if
        end do


      case (cExp)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = exp(Dstack(iblock, istackPtr))
        end do


      case (cFloor)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = floor(Dstack(iblock, istackPtr))
        end do


      case (cIf)
        ! IF-THEN-ELSE cannot be vectorised which should be noted during
        ! bytecode compilation. If we reach this point, then something
        ! went wrong before.
        call output_line('IF-THEN-ELSE cannot be vectorised!',&
            OU_CLASS_ERROR, OU_MODE_STD,'evalFunctionBlock')
        call sys_halt()


      case (cLog)
        do iblock = 1, iblockSize
          if (Dstack(iblock, istackPtr) .le. 0._DP) then
            EvalErrType = 3
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = log(Dstack(iblock, istackPtr))
          end if
        end do


      case (cLog10)
        do iblock = 1, iblockSize
          if (Dstack(iblock, istackPtr) .le. 0._DP) then
            EvalErrType = 3
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = log10(Dstack(iblock, istackPtr))
          end if
        end do


      case (cMax)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = max(Dstack(iblock, istackPtr-1),&
                                             Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cMin)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = min(Dstack(iblock, istackPtr-1),&
                                             Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1

      case (cRrand)
        call random_seed (size=iaux)
        allocate (p_Irandom1(iaux))
        allocate (p_Irandom2(iaux))
        call random_seed (get=p_Irandom1)

        do iblock = 1, iblockSize
          p_Irandom2(:) = 0
          p_Irandom2(1) = int(Dstack(iblock, istackPtr-1))
          call random_seed (put=p_Irandom2)
          daux = 0.0_DP
          do iaux=1,max(1,int(Dstack(iblock, istackPtr)))
            call random_number (daux)
          end do
          Dstack(iblock, istackPtr-1) = daux
        end do

        call random_seed (put=p_Irandom1)
        deallocate(p_Irandom1,p_Irandom2)

        istackPtr = istackPtr-1


      case (cSec)
        do iblock = 1, iblockSize
          daux = cos(Dstack(iblock, istackPtr))
          if (daux .eq. 0._DP) then
            EvalErrType = 1
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = 1._DP/daux
          end if
        end do


      case (cSign)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = sign(1._DP,Dstack(iblock, istackPtr))
        end do


      case (cSin)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = sin(Dstack(iblock, istackPtr))
        end do


      case(cSinh)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = sinh(Dstack(iblock, istackPtr))
        end do


      case(cSqrt)
        do iblock = 1, iblockSize
          if (Dstack(iblock, istackPtr) .lt. 0._DP) then
            EvalErrType = 3
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr) = sqrt(Dstack(iblock, istackPtr))
          end if
        end do


      case (cTan)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = tan(Dstack(iblock, istackPtr))
        end do


      case (cTanh)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = tanh(Dstack(iblock, istackPtr))
        end do


        !------------------------------------------------------------
        ! Misc
        !------------------------------------------------------------
      case (cImmed)
        istackPtr = istackPtr+1
        do iblock = 1, iblockSize
          Dstack(iblock, istackPtr) = rcomp%Dimmed(idataPtr)
        end do
        idataPtr = idataPtr+1


      case (cJump)
        idataPtr = rcomp%IbyteCode(iinstPtr+2)
        iinstPtr = rcomp%IbyteCode(iinstPtr+1)


        !------------------------------------------------------------
        ! Operators
        !------------------------------------------------------------
      case (cNeg)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = -Dstack(iblock, istackPtr)
        end do


      case (cAdd)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = Dstack(iblock, istackPtr-1)+&
                                         Dstack(iblock, istackPtr)
        end do
        istackPtr = istackPtr-1


      case (cSub)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = Dstack(iblock, istackPtr-1)-&
                                         Dstack(iblock, istackPtr)
        end do
        istackPtr = istackPtr-1


      case (cMul)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = Dstack(iblock, istackPtr-1)*&
                                         Dstack(iblock, istackPtr)
        end do
        istackPtr = istackPtr-1


      case (cDiv)
        do iblock = 1, iblockSize
          if (Dstack(iblock, istackPtr) .eq. 0._DP) then
            EvalErrType = 1
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr-1) = Dstack(iblock, istackPtr-1)/&
                                          Dstack(iblock, istackPtr)
          end if
        end do
        istackPtr = istackPtr-1


      case (cMod)
        do iblock = 1, iblockSize
          if (Dstack(iblock, istackPtr) .eq. 0._DP) then
            EvalErrType = 1
            Dresult(iblock) = 0._DP
          else
            Dstack(iblock, istackPtr-1) = mod(Dstack(iblock, istackPtr-1),&
                                              Dstack(iblock, istackPtr))
          end if
        end do
        istackPtr = istackPtr-1


      case (cPow)
        do  iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = Dstack(iblock, istackPtr-1)**Dstack(iblock, istackPtr)
        end do
        istackPtr = istackPtr-1


      case (cEqual)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(Dstack(iblock, istackPtr-1) .eq.&
                                                    Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cNEqual)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(Dstack(iblock, istackPtr-1) .ne.&
                                                    Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cLess)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(Dstack(iblock, istackPtr-1) .lt.&
                                                    Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cLessOrEq)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(Dstack(iblock, istackPtr-1) .le.&
                                                    Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cGreater)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(Dstack(iblock, istackPtr-1) .gt.&
                                                    Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cGreaterOrEq)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(Dstack(iblock, istackPtr-1) .ge.&
                                                    Dstack(iblock, istackPtr))
        end do
        istackPtr = istackPtr-1


      case (cAnd)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(DbleToLogc(Dstack(iblock, istackPtr-1)) .and. &
                                                    DbleToLogc(Dstack(iblock, istackPtr)) )
        end do
        istackPtr = istackPtr-1


      case (cOr)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr-1) = LogcToDble(DbleToLogc(Dstack(iblock, istackPtr-1)) .or. &
                                                    DbleToLogc(Dstack(iblock, istackPtr)) )
        end do
        istackPtr = istackPtr-1


      case (cNot)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = LogcToDble( .not. DbleToLogc(Dstack(iblock, istackPtr)) )
        end do


        !------------------------------------------------------------
        ! Degrees-radians conversion
        !------------------------------------------------------------
      case (cDeg)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = RadToDeg(Dstack(iblock, istackPtr))
        end do


      case (cRad)
        do iblock = 1, iblockSize
           Dstack(iblock, istackPtr) = DegToRad(Dstack(iblock, istackPtr))
        end do


      case default
        istackPtr  = istackPtr+1
        iVariable = rcomp%IbyteCode(iinstPtr)-VarBegin+1

        ! Do we have to process one of the scalar variables of one of the block variables
        if (iVariable .ge. istartValueScalar) then
          do iblock = 1, iblockSize
             Dstack(iblock, istackPtr) = DvalueScalar(iVariable-istartValueScalar+1)
          end do
        else
          if (idim .eq. 1) then
            do iblock = 1, iblockSize
               Dstack(iblock, istackPtr) = DvalueBlock(iblock, iVariable)
            end do
          else
            do iblock = 1, iblockSize
               Dstack(iblock, istackPtr) = DvalueBlock(iVariable, iblock)
            end do
          end if
        end if
      end select
    end do

    EvalErrType = 0
    do iblock = 1, iblockSize
       Dresult(iblock) = Dstack(iblock, istackPtr)
    end do

  end subroutine evalFunctionBlock

end module fparser
