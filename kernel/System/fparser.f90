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

MODULE fparser
  USE fsystem
  USE storage
  USE stack
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_fparser
  PUBLIC :: fparser_init
  PUBLIC :: fparser_done
  PUBLIC :: fparser_defineConstant
  PUBLIC :: fparser_defineExpression
  PUBLIC :: fparser_create
  PUBLIC :: fparser_release
  PUBLIC :: fparser_parseFunction
  PUBLIC :: fparser_evalFunction
  PUBLIC :: fparser_ErrorMsg
  PUBLIC :: fparser_PrintByteCode
  
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  
  INTERFACE fparser_evalFunction
    MODULE PROCEDURE fparser_evalFunctionScalar
    MODULE PROCEDURE fparser_evalFunctionBlock
  END INTERFACE

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<constants>

!<constantblock description="Global constants for parser">
  ! Maximum stack size
  ! The stack is used for all instances of function parsers. During initialization
  ! a smaller size can be specified. This is only the maximum size of the stack.
  INTEGER, PARAMETER :: FPAR_MAXSTACKSIZE   = 1024*1024 ! this is 8 MByte

  ! Lenght of string
  INTEGER, PARAMETER :: FPAR_STRLEN         = 2048

  ! Maximum number of predefined/user-defined constants
  INTEGER, PARAMETER :: FPAR_MAXCONSTS      = 128

  ! Maximum number of predefined/user-defined expressions
  INTEGER, PARAMETER :: FPAR_MAXEXPRESSIONS = 128

  ! Length of constant name
  INTEGER, PARAMETER :: FPAR_CONSTLEN       = 12

  ! Length of expression name
  INTEGER, PARAMETER :: FPAR_EXPRLEN        = 12
!</constantblock>


!<constantblock description="data type for parser bytecode">
  ! Data type of bytecode
  INTEGER, PARAMETER :: is = SELECTED_INT_KIND(1)
!</constantblock>


!<constantblock description="keywords for parser">
  INTEGER(is), PARAMETER :: cImmed       =  1, &
                            cJump        =  2, &
                            cNeg         =  3, &
                            cDeg         =  4, &
                            cRad         =  5, &
                            cAdd         =  6, & ! <-- first dyadic operator: A .OP. B
                            cSub         =  7, &
                            cMul         =  8, &
                            cDiv         =  9, &
                            cMod         = 10
  INTEGER(is), PARAMETER :: cPow         = 11, &
                            cNEqual      = 12, & ! NOTE: != must be prior to =
                            cEqual       = 13, &
                            cLessOrEq    = 14, & ! NOTE: <= must be prior to < 
                            cLess        = 15, &
                            cGreaterOrEq = 16, & ! NOTE: >= must be prior to > 
                            cGreater     = 17, &
                            cNot         = 18, &
                            cAnd         = 19, &
                            cOr          = 20    ! --> last dyadic operator: A.OP. B
  INTEGER(is), PARAMETER :: cIf          = 21, & ! <-- if-then-else
                            cMin         = 22, & ! <-- first dyadic operator: .OP.(A,B)
                            cMax         = 23, &
                            cAtan2       = 24, & ! --> last dyadic operator: .OP.(A,B)
                            cAbs         = 25, & ! <-- monadic operator: .OP.(A)
                            cAnint       = 26, &
                            cAint        = 27, &
                            cExp         = 28, &
                            cLog10       = 29, &
                            cLog         = 30
  INTEGER(is), PARAMETER :: cSqrt        = 31, &
                            cSinh        = 32, &
                            cCosh        = 33, &
                            cTanh        = 34, &
                            cSin         = 35, &
                            cCos         = 36, &
                            cTan         = 37, & 
                            cCot         = 38, &
                            cAsinh       = 39, &
                            cAsin        = 40
  INTEGER(is), PARAMETER :: cAcosh       = 41, &
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
  CHARACTER (LEN=2), DIMENSION(cAdd:cOr), PARAMETER :: Ops = (/ '+ ', &
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
  CHARACTER (LEN=5), DIMENSION(cIf:cSec), PARAMETER :: Funcs = (/ 'if   ', &
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
  CHARACTER(LEN=FPAR_CONSTLEN), DIMENSION(3) :: PredefinedConsts = (/ 'pi        ', &
                                                                      'exp       ', &
                                                                      'infty     ' /)
!</constantblock>


!<constantblock description="predefined constant values for parser">
  REAL(DP), DIMENSION(3), PARAMETER :: PredefinedConstvals = (/&
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
  INTEGER, SAVE                                                    :: h_Stack = ST_NOHANDLE
  
  ! Global pointer to stack memory
  REAL(DP), DIMENSION(:), POINTER, SAVE                            :: p_Stack => NULL()

  ! Global number of predefined/user-defined constants
  INTEGER, SAVE                                                    :: nConsts = 0 

  ! Global constant names for parser
  CHARACTER(LEN=FPAR_CONSTLEN), DIMENSION(FPAR_MAXCONSTS), SAVE    :: Consts  = '     '

  ! Global constant values for parser
  REAL(DP), DIMENSION(FPAR_MAXCONSTS), SAVE                        :: ConstVals = 0

  ! Global number of predefined/user-defined expressions
  INTEGER, SAVE                                                    :: nExpressions = 0

  ! Global expression name for parser
  CHARACTER(LEN=FPAR_EXPRLEN), DIMENSION(FPAR_MAXCONSTS), SAVE     :: Expressions = ''

  ! Global expression string for parser
  CHARACTER(LEN=FPAR_STRLEN), DIMENSION(FPAR_MAXEXPRESSIONS), SAVE :: ExpressionStrings
!</globals>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<types>

!<typeblock>
  
  ! Type block for storing all information of the function parser
  TYPE t_fparser
    PRIVATE
    
    ! Array of function parser components. Each component is used to
    ! handle one function string at a time
    TYPE(t_fparserComponent), DIMENSION(:), POINTER :: Comp => NULL()
    
    ! Number of parser components
    INTEGER :: nComp                                        = 0
  END TYPE t_fparser

!</typeblock>

!<typeblock>

  ! Type block for storing the bytecode of the function parser for one component
  TYPE t_fparserComponent
    PRIVATE

    ! Size of bytecode
    INTEGER :: ByteCodeSize                                 = 0

    ! Size of immediates
    INTEGER :: ImmedSize                                    = 0

    ! Stack size
    INTEGER :: StackSize                                    = 0

    ! Stack pointer
    INTEGER :: StackPtr                                     = 0

    ! Use degree conversion DEG <-> RAD for some functions
    LOGICAL :: useDegreeConversion                          = .FALSE.

    ! Is vectorizable
    LOGICAL :: isVectorizable                               = .TRUE.
    
    ! Bytecode
    INTEGER(is), DIMENSION(:), POINTER :: ByteCode          => NULL()
    
    ! Immediates
    REAL(DP), DIMENSION(:), POINTER :: Immed                => NULL()
  END TYPE t_fparserComponent
!</typeblock>

!</types>

CONTAINS

  ! *****************************************************************************
  
!<subroutine>

  SUBROUTINE fparser_init (iStackSize)

!<description>
    ! Initialize function parser
!</description>

!<input>
    ! OPTIONAL: initial size of the stack memory
    INTEGER, INTENT(IN), OPTIONAL :: iStackSize
!</input>
!</subroutine>

    INTEGER :: iSize,iConst

    IF (PRESENT(iStackSize)) THEN
      iSize=MIN(iStackSize,FPAR_MAXSTACKSIZE)
    ELSE
      iSize=FPAR_MAXSTACKSIZE
    END IF
    
    ! Allocate memory for global stack and set pointer
    CALL storage_new('fparser_init','p_Stack',iSize,ST_DOUBLE,h_Stack,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_double(h_Stack,p_Stack)
    
    ! Initialize predefined constants
    DO iConst=LBOUND(PredefinedConsts,1),UBOUND(PredefinedConsts,1)
      CALL fparser_defineConstant(PredefinedConsts(iConst),PredefinedConstvals(iConst))
    END DO
  END SUBROUTINE fparser_init

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_done

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
    IF (h_Stack  /= ST_NOHANDLE) CALL storage_free(h_Stack)
    NULLIFY(p_Stack)
  END SUBROUTINE fparser_done

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_defineConstant(constant,constantValue)

!<description>
    ! Define a new constant for the function parser.
    ! This subroutine checks if the given constant is already defined.
!</description>

!<input>
    ! Name of the constant
    CHARACTER(LEN=FPAR_CONSTLEN), INTENT(IN) :: constant

    ! Value of the constant
    REAL(DP), INTENT(IN) :: constantValue
!</input>
!</subroutine>

    ! local variables
    CHARACTER(LEN=LEN(constant)) :: const
    INTEGER :: iConst

    ! Check if there is enough space
    IF (nConsts < FPAR_MAXCONSTS) THEN
      nConsts = nConsts+1
    ELSE
      PRINT *, "*** Parser error: no space left for definition of constant"
      CALL sys_halt()
    END IF

    ! Prepare constant
    CALL sys_tolower(constant,const)

    ! Check if constant is already defined
    DO iConst=1,nConsts-1
      IF (Consts(iConst) == const) THEN
        ! If it is already defined, then it must not have a different value
        IF(ConstVals(iConst) /= constantValue) THEN
          PRINT *, "*** Parser error: constant is already defined with different value"
          CALL sys_halt()
        ELSE
          nConsts = nConsts-1
          RETURN
        END IF
      END IF
    END DO

    ! Apply constant value and constant name
    Consts(nConsts)    = const
    ConstVals(nConsts) = constantValue
  END SUBROUTINE fparser_defineConstant

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_defineExpression(expression,expressionString)

!<description>
    ! Define a new expression for the function parser.
    ! This subroutine checks if the given expression is already defined.
!</description>

!<input>
    ! Name of the expression
    CHARACTER(LEN=FPAR_EXPRLEN), INTENT(IN) :: expression

    ! String of the expression
    CHARACTER(LEN=*), INTENT(IN) :: expressionString
!</input>
!</subroutine>

    ! local variables
    CHARACTER(LEN=LEN(expression))       :: expr
    CHARACTER(LEN=LEN(expressionString)) :: str
    INTEGER :: iExpression

    ! Check if there is enough space
    IF (nExpressions < FPAR_MAXEXPRESSIONS) THEN
      nExpressions = nExpressions+1
    ELSE
      PRINT *, "*** Parser error: no space left for definition of expression"
      CALL sys_halt()
    END IF

    ! Prepare expression string
    CALL sys_tolower(expression,expr)
    CALL sys_tolower(expressionString,str)

    ! Replace human readable function names by 1-Char. format
    CALL Replace ('**','^ ',str)
    CALL Replace ('[','(',  str)
    CALL Replace (']',')',  str)
    CALL Replace ('{','(',  str)
    CALL Replace ('}',')',  str)
    
    ! Condense function string
    CALL RemoveSpaces (str)

    ! Check if expressions is already defined
    DO iExpression=1,nExpressions-1
      IF (Expressions(iExpression) == expr) THEN
        ! If it is already defined, then it must not have a different value
        IF(ExpressionStrings(iExpression) /= str) THEN
          PRINT *, "*** Parser error: expression is already defined with different string"
          CALL sys_halt()
        ELSE
          nExpressions = nExpressions-1
          RETURN
        END IF
      END IF
    END DO

    ! Apply expressions string and expression name
    Expressions(nExpressions)       = expr
    ExpressionStrings(nExpressions) = str
  END SUBROUTINE fparser_defineExpression

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE fparser_create (rparser,nComp)

!<description>
    ! Initialize function parser for nComp functions
!</description>

!<input>
    ! Number of functions
    INTEGER, INTENT(IN) :: nComp
!</input>

!<inputoutput>
    ! Function parser object
    TYPE(t_fparser), INTENT(INOUT) :: rparser
!</inputoutput>
!</subroutine>
    
    ! Set number of components
    rparser%nComp=nComp
    
    IF (ASSOCIATED(rparser%Comp)) DEALLOCATE(rparser%Comp)
    ALLOCATE (rparser%Comp(nComp))
  END SUBROUTINE fparser_create

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_release (rparser)

!<description>
    ! Release function parser and all of its coponents
!</description>

!<inputoutput>
    ! Function parser
    TYPE(t_fparser), INTENT(INOUT) :: rparser
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER :: iComp
    
    ! Check that pointer is associated and return otherwise
    IF (.NOT.ASSOCIATED(rparser%Comp)) RETURN

    ! Loop over all components
    DO iComp=1,rparser%nComp
      IF (ASSOCIATED(rparser%Comp(iComp)%ByteCode)) DEALLOCATE(rparser%Comp(iComp)%ByteCode)
      IF (ASSOCIATED(rparser%Comp(iComp)%Immed))    DEALLOCATE(rparser%Comp(iComp)%Immed)
    END DO
    DEALLOCATE(rparser%Comp)
    rparser%nComp=0
  END SUBROUTINE fparser_release

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE fparser_parseFunction (rparser, iComp, FuncStr, Var, useDegrees)

!<description>
    ! Parse ith function string FuncStr and compile it into bytecode
!</description>
    
!<input>
    ! Function identifier
    INTEGER, INTENT(IN) :: iComp

    ! Function string
    CHARACTER (LEN=*), INTENT(IN) :: FuncStr

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: Var

    ! OPTIONAL
    LOGICAL, INTENT(IN), OPTIONAL :: useDegrees
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_fparser), INTENT(INOUT) :: rparser
!</inputoutput>
!</subroutine>

    ! local variables
    CHARACTER (LEN=LEN(FuncStr)) :: Func
    
    IF (iComp < 1 .OR. iComp > rparser%nComp) THEN
      WRITE(*,*) '*** Parser error: Function number ',iComp,' out of range'
      CALL sys_halt()
    END IF
    
    ! Local copy of function string
    Func = FuncStr

    ! Replace human readable function names by 1-Char. format
    CALL Replace ('**','^ ',Func)
    CALL Replace ('[','(',Func)
    CALL Replace (']',')',Func)
    CALL Replace ('{','(',Func)
    CALL Replace ('}',')',Func)

    ! Condense function string
    CALL RemoveSpaces (Func)

    ! Check for valid syntax; this prevents the bytecode compiler
    ! from running into endless loops or other problems
    CALL CheckSyntax (Func,Var)
    
    ! If syntax is correct, then compile into bytecode
    IF (PRESENT(useDegrees))&
        rparser%Comp(iComp)%useDegreeConversion = useDegrees

    CALL Compile (rparser%Comp(iComp),Func,Var)
  END SUBROUTINE fparser_parseFunction

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE fparser_evalFunctionScalar (rparser, iComp, Val, Res)

!<description>
    ! Evaluate bytecode of component iComp for the values passed in array 
    ! Val(:). Note that this function is a wrapper for the working routine
    ! evalFunctionScalar. It is used to adjust the dimensions of the global
    ! stack memory if required.
!</description>

!<input>
    ! Function identifier
    INTEGER, INTENT(IN) :: iComp

    ! Variable values
    REAL(DP), DIMENSION(:), INTENT(IN) :: Val
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_fparser),  INTENT(IN) :: rparser
!</inputoutput>

!<output>
    ! Evaluated function
    REAL(DP), INTENT(OUT)  :: Res
!</output>
!</subroutine>

    ! local variables
    INTEGER :: EvalErrType

    IF (h_Stack .EQ. ST_NOHANDLE) THEN
      PRINT *, "*** Parser error: Parser not initialised!"
      CALL sys_halt()
    END IF
    
    ! Check if memory of global stack is sufficient
    IF (SIZE(p_Stack) < rparser%Comp(iComp)%StackSize+1) THEN
      IF (rparser%Comp(iComp)%StackSize+1 < FPAR_MAXSTACKSIZE) THEN
        CALL storage_realloc('fparser_evalFunctionScalar',&
            rparser%Comp(iComp)%StackSize+1,h_Stack,ST_NEWBLOCK_NOINIT,.FALSE.)
        CALL storage_getbase_double(h_Stack,p_Stack)
      ELSE
        PRINT *, "*** Parser error: stack size exceeds memory!"
        CALL sys_halt()
      END IF
    END IF
    
    ! Invoke working routine
    CALL evalFunctionScalar(p_Stack,rparser%Comp(iComp),Val,EvalErrType,Res)

    ! Check if evaluation was successful
    IF (EvalErrType .NE. 0) THEN
      PRINT *, "*** Parser error: An error occured during function evaluation!"
      CALL sys_halt()
    END IF
  END SUBROUTINE fparser_evalFunctionScalar

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_evalFunctionBlock (rparser, iComp, iDim, ValBlock, Res, ValScalar)

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
    ! Function identifier
    INTEGER, INTENT(IN) :: iComp

    ! Orientation of the stored values
    ! iDim =1 : ValBlock is organized as (x1:xN),(y1:yN),...
    ! iDim =2 : ValBlock is organized as (x1,y1),(x2,y2),...,(xN,yN)
    INTEGER, INTENT(IN) :: iDim

    ! Variable values (must have the same dimension as Res)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: ValBlock
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_fparser),  INTENT(INOUT) :: rparser

    ! Variable values. This is a vector of scalar variables
    ! which is the same for all components of Res, e.g. the time variable.
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: ValScalar
!</inputoutput>

!<output>
    ! Evaluated function
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Res
!</output>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ValTemp
    INTEGER :: iVal,jVal,nVal,iMemory,iBlockSize,iBlock,sizeValBlock,sizeValScalar
    INTEGER :: EvalErrType

    ! Get total number of variable sets
    nVal=SIZE(ValBlock,iDim)

    ! Check if the compiled function is vectorizable
    IF (rparser%Comp(iComp)%isVectorizable) THEN
      ! ...ok, vectorization of the bytecode is admissible.

      ! Guess required memory
      iMemory=(rparser%Comp(iComp)%StackSize+1)*nVal
      iMemory=MIN(iMemory,FPAR_MAXSTACKSIZE)
      
      ! Check if memory of global stack is sufficient for complete set of variables
      IF (SIZE(p_Stack) < iMemory) THEN
        IF (rparser%Comp(iComp)%StackSize+1 < FPAR_MAXSTACKSIZE) THEN
          CALL storage_realloc('fparser_evalFunctionBlock',&
              iMemory,h_Stack,ST_NEWBLOCK_NOINIT,.FALSE.)
          CALL storage_getbase_double(h_Stack,p_Stack)
        ELSE
          PRINT *, "*** Parser error: stack size exceeds memory!"
          CALL sys_halt()
        END IF
      END IF
      
      ! Compute size of blocks that can be vectorized
      iBlockSize = FLOOR(REAL(iMemory)/REAL(rparser%Comp(iComp)%StackSize+1))
      
      ! What is the organization of ValBlock(:,:)
      IF (iDim == 1) THEN
        DO iVal=1,nVal,iBlockSize
          
          ! Initialization
          jVal     = MIN(iVal+iBlockSize,nVal)
          iBlock   = jVal-iVal+1
          
          ! Invoke working routine
          CALL evalFunctionBlock(iBlock,p_Stack,rparser%Comp(iComp),&
              ValBlock(iVal:jVal,:),iDim,EvalErrType,Res(iVal:jVal),ValScalar)
        END DO
      ELSE
        DO iVal=1,nVal,iBlockSize
          
          ! Initialization
          jVal     = MIN(iVal+iBlockSize,nVal)
          iBlock   = jVal-iVal+1
          
          ! Invoke working routine
          CALL evalFunctionBlock(iBlock,p_Stack,rparser%Comp(iComp),&
              ValBlock(:,iVal:jVal),iDim,EvalErrType,Res(iVal:jVal),ValScalar)
        END DO
      END IF
      
    ELSE   ! The compiled function cannot be vectorized

      ! Check if memory of global stack is sufficient for one set of variables
      IF (SIZE(p_Stack) < rparser%Comp(iComp)%StackSize+1) THEN
        IF (rparser%Comp(iComp)%StackSize+1 < FPAR_MAXSTACKSIZE) THEN
          CALL storage_realloc('fparser_evalFunctionBlock',&
              rparser%Comp(iComp)%StackSize+1,h_Stack,ST_NEWBLOCK_NOINIT,.FALSE.)
          CALL storage_getbase_double(h_Stack,p_Stack)
        ELSE
          PRINT *, "*** Parser error: stack size exceeds memory!"
          CALL sys_halt()
        END IF
      END IF

      ! The compiled bytecode cannot be vectorized. Hence, evaluate the function
      ! separately for each set of variables. Here, the organization of the array
      ! VAL(:,:) is important. Moreover, if the optional parameter ValScalar is
      ! given, then we have to combine those variables from ValBlock and ValScalar.

      IF (PRESENT(ValScalar)) THEN

        ! Allocate auxiliary array
        sizeValBlock = SIZE(ValBlock,3-iDim)
        sizeValScalar= SIZE(ValScalar)
        ALLOCATE(ValTemp(sizeValBlock+sizeValScalar))

        IF (iDim == 1) THEN
          DO iVal=1,nVal
            
            ValTemp(1:sizeValBlock) = ValBlock(iVal,:)
            ValTemp(sizeValBlock+1:)= ValScalar

            ! Invoke working routine
            CALL evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValTemp,&
                EvalErrType,Res(iVal))
          END DO
        ELSE
          DO iVal=1,nVal
            
            ValTemp(1:sizeValBlock) = ValBlock(:,iVal)
            ValTemp(sizeValBlock+1:)= ValScalar
            
            ! Invoke working routine
            CALL evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValTemp,&
                EvalErrType,Res(iVal))
          END DO
        END IF

        ! Deallocate auxiliary array
        DEALLOCATE(ValTemp)

      ELSE
        
        IF (iDim == 1) THEN
          DO iVal=1,nVal
            
            ! Invoke working routine
            CALL evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValBlock(iVal,:),&
                EvalErrType,Res(iVal))
          END DO
        ELSE
          DO iVal=1,nVal
            
            ! Invoke working routine
            CALL evalFunctionScalar(p_Stack,rparser%Comp(iComp),ValBlock(:,iVal),&
                EvalErrType,Res(iVal))
          END DO
        END IF

      END IF
    END IF

    ! Check if evaluation was successful
    IF (EvalErrType .NE. 0) THEN
      PRINT *, "*** Parser error: An error occured during function evaluation!"
      CALL sys_halt()
    END IF
  END SUBROUTINE fparser_evalFunctionBlock

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE evalFunctionScalar (Stack, Comp, Val, EvalErrType, res)

!<description>
    ! Evaluate bytecode for the values passed in array Val(:). Note, this subroutine
    ! requires some working memory Stack(*) which is not checked. So be warned,
    ! not to call this routine directly ;-)
!</description>

!<input>
    ! Component of function parser
    TYPE(t_fparserComponent), INTENT(IN) :: Comp

    ! Variable values
    REAL(DP), DIMENSION(:), INTENT(IN) :: Val
!</input>

!<inputoutput>
    ! Stack memory
    REAL(DP), DIMENSION(*), INTENT(INOUT) :: Stack
!</inputoutput>

!<output>
    ! Error code for function evaluation
    INTEGER, INTENT(OUT) :: EvalErrType

    ! Evaluated function
    REAL(DP), INTENT(OUT) :: res
!</output>
!</subroutine>

    ! local parameters
    REAL(DP), PARAMETER :: Zero = 0._DP
    
    ! local variables
    INTEGER  :: InstPtr,StackPtr,DataPtr
    INTEGER  :: jumpAddr,immedAddr
    REAL(DP) :: daux
    
    ! Initialization
    DataPtr  = 1
    StackPtr = 0
    InstPtr  = 0
    
    ! Repeat until complete bytecode has been processed
    DO WHILE(InstPtr < Comp%ByteCodeSize)
      InstPtr=InstPtr+1

      ! What kind of bytecode are we?
      SELECT CASE (Comp%ByteCode(InstPtr))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      CASE (cAbs)
        Stack(StackPtr)=ABS(Stack(StackPtr))
        
      CASE (cAcos)
        IF ((Stack(StackPtr) < -1._DP) .OR. &
            (Stack(StackPtr) >  1._DP)) THEN
          EvalErrType = 4
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr)=ACOS(Stack(StackPtr))
        
      CASE (cAsin)
        IF ((Stack(StackPtr) < -1._DP) .OR. &
            (Stack(StackPtr) >  1._DP)) THEN
          EvalErrType = 4
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr)=ASIN(Stack(StackPtr))
        
      CASE (cAtan)
        Stack(StackPtr)=ATAN(Stack(StackPtr))

      CASE (cAtan2)
        Stack(StackPtr-1)=ATAN2(Stack(StackPtr -1),Stack(StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cAcosh)
        daux=Stack(StackPtr)+SQRT(Stack(StackPtr)**2-1)
        IF (daux <= 0._DP) THEN
          EvalErrType = 5
          res         = zero
          RETURN
        END IF
        Stack(StackPtr)=LOG(daux)
        
      CASE (cAnint)
        Stack(StackPtr)=ANINT(Stack(StackPtr))

      CASE (cAint)
        Stack(StackPtr)=AINT(Stack(StackPtr))

      CASE (cAsinh)
        daux=Stack(StackPtr)+SQRT(Stack(StackPtr)**2+1)
        IF (daux <= 0._DP) THEN
          EvalErrType = 5
          res         = zero
          RETURN
        END IF
        Stack(StackPtr)=LOG(daux)
        
      CASE (cAtanh)
        IF (Stack(StackPtr) == -1._DP) THEN
          EvalErrType = 6
          res         = zero
          RETURN
        END IF
        daux=(1+Stack(StackPtr))/(1-Stack(StackPtr))
        IF (daux <= 0._DP) THEN
          EvalErrType = 3
          res         = zero
          RETURN
        END IF
        Stack(StackPtr)=LOG(daux)/2._DP
        
      CASE (cCeil)
        Stack(StackPtr)=CEILING(Stack(StackPtr))
        
      CASE (cCos)
        Stack(StackPtr)=COS(Stack(StackPtr)) 
        
      CASE (cCosh)
        Stack(StackPtr)=COSH(Stack(StackPtr))

      CASE (cCot)
        daux=TAN(Stack(StackPtr))
        IF (daux == 0._DP) THEN 
          EvalErrType = 1
          res         = zero
          RETURN
        END IF
        Stack(StackPtr)=1._DP/daux

      CASE (cCsc)
        daux=SIN(Stack(StackPtr))
        IF (daux == 0._DP) THEN
          EvalErrType = 1
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr)=1._DP/daux
        
      CASE (cExp)
        Stack(StackPtr)=EXP(Stack(StackPtr))

      CASE (cFloor)
        Stack(StackPtr)=FLOOR(Stack(StackPtr))

      CASE (cIf)
        InstPtr=InstPtr+1; jumpAddr  = Comp%ByteCode(InstPtr)
        InstPtr=InstPtr+1; immedAddr = Comp%ByteCode(InstPtr)
        IF (.NOT.DbleToLogc(Stack(StackPtr))) THEN
          InstPtr = jumpAddr
          DataPtr = immedAddr
        END IF
        StackPtr=StackPtr-1
        
      CASE (cLog)
        IF (Stack(StackPtr) <= 0._DP) THEN
          EvalErrType = 3
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr)=LOG(Stack(StackPtr)) 
        
      CASE (cLog10)
        IF (Stack(StackPtr) <= 0._DP) THEN
          EvalErrType = 3
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr)=LOG10(Stack(StackPtr))
        
      CASE (cMax)
        Stack(StackPtr-1)=MAX(Stack(StackPtr-1),Stack(StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cMin)
        Stack(StackPtr-1)=MIN(Stack(StackPtr-1),Stack(StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cSec)
        daux=COS(Stack(StackPtr))
        IF (daux == 0._DP) THEN
          EvalErrType = 1
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr)=1._DP/daux
        
      CASE (cSin)
        Stack(StackPtr)=SIN(Stack(StackPtr))
        
      CASE(cSinh)
        Stack(StackPtr)=SINH(Stack(StackPtr))
      
      CASE(cSqrt)
        IF (Stack(StackPtr) < 0._DP) THEN
          EvalErrType = 3
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr)=SQRT(Stack(StackPtr))

      CASE (cTan)
        Stack(StackPtr)=TAN(Stack(StackPtr))
        
      CASE (cTanh)
        Stack(StackPtr)=TANH(Stack(StackPtr))
        
        !------------------------------------------------------------
        ! Misc
        !------------------------------------------------------------
      CASE (cImmed)
        StackPtr        = StackPtr+1
        Stack(StackPtr) = Comp%Immed(DataPtr)
        DataPtr         = DataPtr+1
        
      CASE (cJump)
        DataPtr=Comp%ByteCode(InstPtr+2)
        InstPtr=Comp%ByteCode(InstPtr+1)

        !------------------------------------------------------------
        ! Operators
        !------------------------------------------------------------
      CASE (cNeg)
        Stack(StackPtr)=-Stack(StackPtr)
        
      CASE (cAdd)
        Stack(StackPtr-1)=Stack(StackPtr-1)+Stack(StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cSub)
        Stack(StackPtr-1)=Stack(StackPtr-1)-Stack(StackPtr)
        StackPtr=StackPtr-1

      CASE (cMul)
        Stack(StackPtr-1)=Stack(StackPtr-1)*Stack(StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cDiv)
        IF (Stack(StackPtr) == 0._DP) THEN
          EvalErrType = 1
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr-1)=Stack(StackPtr-1)/Stack(StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cMod)
        IF (Stack(StackPtr) == 0._DP) THEN
          EvalErrType = 1
          res         = zero
          RETURN
        ENDIF
        Stack(StackPtr-1)=MOD(Stack(StackPtr-1),Stack(StackPtr))
        StackPtr=StackPtr-1

      CASE (cPow)
        Stack(StackPtr-1)=Stack(StackPtr-1)**Stack(StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cEqual)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) == Stack(StackPtr))
        StackPtr=StackPtr-1

      CASE (cNEqual)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) /= Stack(StackPtr))
        StackPtr=StackPtr-1

      CASE (cLess)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) < Stack(StackPtr))
        StackPtr=StackPtr-1

      CASE (cLessOrEq)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) <= Stack(StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cGreater)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) > Stack(StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cGreaterOrEq)
        Stack(StackPtr-1)=LogcToDble(Stack(StackPtr-1) >= Stack(StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cAnd)
        Stack(StackPtr-1)=LogcToDble(DbleToLogc( Stack(StackPtr-1)) .AND. &
            DbleToLogc(Stack(StackPtr)) )
        StackPtr=StackPtr-1

      CASE (cOr)
        Stack(StackPtr-1)=LogcToDble(DbleToLogc( Stack(StackPtr-1)) .OR. &
            DbleToLogc(Stack(StackPtr)) )
        StackPtr=StackPtr-1

      CASE (cNot)
        Stack(StackPtr)=LogcToDble( .NOT. DbleToLogc(Stack(StackPtr)) )
        
        !------------------------------------------------------------
        ! Degrees-radians conversion
        !------------------------------------------------------------
      CASE (cDeg)
        Stack(StackPtr)=RadToDeg(Stack(StackPtr))
        
      CASE (cRad)
        Stack(StackPtr)=DegToRad(Stack(StackPtr))
        
      CASE DEFAULT
        StackPtr=StackPtr+1
        Stack(StackPtr)=Val(Comp%ByteCode(InstPtr)-VarBegin+1)
      END SELECT
    END DO
    
    EvalErrType = 0
    Res = Stack(StackPtr)
  END SUBROUTINE evalFunctionScalar

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE evalFunctionBlock (BlockSize, Stack, Comp, ValBlock, iDim, EvalErrType, Res, ValScalar)

!<description>
    ! Evaluate bytecode for an array of values passed in ValBlock(:,:).
    ! Note, this subroutine requires some working memory Stack(iBlock,*) which is
    ! not checked. So be warned, not to call this function directly ;-)
!</description>

!<input>
    ! Component of function parser
    TYPE(t_fparserComponent), INTENT(IN) :: Comp
    
    ! Variable values
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: ValBlock

    ! Size of the vector block
    INTEGER, INTENT(IN) :: BlockSize

    ! Orientation of the stored values
    ! iDim =1 : ValBlock is organized as (x1:xN),(y1:yN),...
    ! iDim =2 : ValBlock is organized as (x1,y1),(x2,y2),...,(xN,yN)
    INTEGER, INTENT(IN) :: iDim

    ! Vector of scalar variable values
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: ValScalar
!</input>

!<inputoutput>
    ! Stack memory
    REAL(DP), DIMENSION(BlockSize,*), INTENT(INOUT) :: Stack
!</inputoutput>

!<output>
    ! Error code for function evaluation
    INTEGER, INTENT(OUT) :: EvalErrType

    ! Evaluated function
    REAL(DP), DIMENSION(:) :: Res
!</output>
!</subroutine>

    ! local parameters
    REAL(DP), PARAMETER :: Zero = 0._DP
    
    ! local variables
    INTEGER  :: InstPtr,DataPtr,StackPtr,iBlock,istartValScalar,iVariable
    REAL(DP) :: daux
    
    ! Initialization
    DataPtr  = 1
    StackPtr = 0
    InstPtr  = 0

    ! This is tricky. istartValScalar indicates the number of the first
    ! variable which is passed as scalar. Hence, if the optional parameter
    ! ValScalar is missing, then istartValScalar pointers to SIZE(ValBlock)+1.
    ! Obviously, no variable beyond this value is addressed.
    istartValScalar=SIZE(ValBlock,3-iDim)+1
    
    ! Repeat until complete bytecode has been processed
    DO WHILE(InstPtr < Comp%ByteCodeSize)
      InstPtr=InstPtr+1
      
      SELECT CASE (Comp%ByteCode(InstPtr))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      CASE (cAbs)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=ABS(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        

      CASE (cAcos)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
          IF ((Stack(iBlock,StackPtr) < -1._DP) .OR.&
              (Stack(iBlock,StackPtr) >  1._DP)) THEN
            EvalErrType=4
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr)=ACOS(Stack(iBlock,StackPtr))
          END IF
        END DO
!$omp end parallel do


      CASE (cAsin)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
          IF ((Stack(iBlock,StackPtr) < -1._DP) .OR.&
              (Stack(iBlock,StackPtr) >  1._DP)) THEN
            EvalErrType=4
            Res(iBlock)=zero
          ELSE
            Stack(iBlock,StackPtr)=ASIN(Stack(iBlock,StackPtr))
          END IF
        END DO
!$omp end parallel do

        
      CASE (cAtan)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=ATAN(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do


      CASE (cAtan2)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=ATAN2(Stack(iBlock,StackPtr -1),Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1

        
      CASE (cAcosh)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
          daux=Stack(iBlock,StackPtr)+SQRT(Stack(iBlock,StackPtr)**2-1)
          IF (daux <= 0._DP) THEN
            EvalErrType=5
            Res(iBlock)=zero
          ELSE
            Stack(iBlock,StackPtr)=LOG(daux)
          END IF
        END DO
!$omp end parallel do
        

      CASE (cAnint)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=ANINT(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do


      CASE (cAint)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=AINT(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do


      CASE (cAsinh)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
          daux=Stack(iBlock,StackPtr)+SQRT(Stack(iBlock,StackPtr)**2+1)
          IF (daux <= 0._DP) THEN
            EvalErrType=5
            Res(iBlock)=zero
          ELSE
            Stack(iBlock,StackPtr)=LOG(daux)
          END IF
        END DO
!$omp end parallel do


      CASE (cAtanh)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
          IF (Stack(iBlock,StackPtr) == -1._DP) THEN
            EvalErrType=6
            Res(iBlock)=Zero
          END IF
          daux=(1+Stack(iBlock,StackPtr))/(1-Stack(iBlock,StackPtr))
          IF (daux <= 0._DP) THEN
            EvalErrType=3
            Res(iBlock)=zero
          ELSE
            Stack(iBlock,StackPtr) = LOG(daux)/2._DP
          END IF
        END DO
!$omp end parallel do

        
      CASE (cCeil)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=CEILING(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do


      CASE (cCos)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=COS(Stack(iBlock,StackPtr)) 
        END DO
!$omp end parallel do
        

      CASE (cCosh)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=COSH(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do


      CASE (cCot)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
          daux=TAN(Stack(iBlock,StackPtr))
          IF (daux == 0) THEN 
            EvalErrType=1
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr)=1._DP/daux
          END IF
        END DO
!$omp end parallel do


      CASE (cCsc)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
          daux=SIN(Stack(iBlock,StackPtr))
          IF (daux==0._DP) THEN
            EvalErrType=1
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr)=1._DP/daux
          END IF
        END DO
!$omp end parallel do


      CASE (cExp)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=EXP(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do


      CASE (cFloor)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=FLOOR(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do


      CASE (cIf)
        ! IF-THEN-ELSE cannot be vectorized which should be noted during
        ! bytecode compilation. If we reach this point, then something
        ! went wrong before.
        PRINT *, "*** Parser error: IF-THEN-ELSE cannot be vectorized!"
        CALL sys_halt()
        

      CASE (cLog)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
          IF (Stack(iBlock,StackPtr) <= 0._DP) THEN
            EvalErrType=3
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr)=LOG(Stack(iBlock,StackPtr)) 
          END IF
        END DO
!$omp end parallel do


      CASE (cLog10)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
          IF (Stack(iBlock,StackPtr) <= 0._DP) THEN
            EvalErrType=3
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr)=LOG10(Stack(iBlock,StackPtr))
          END IF
        END DO
!$omp end parallel do


      CASE (cMax)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=MAX(Stack(iBlock,StackPtr-1),Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1
        

      CASE (cMin)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=MIN(Stack(iBlock,StackPtr-1),Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1
        

      CASE (cSec)
!$omp parallel do default(shared) private(iBlock,daux)
        DO iBlock=1,BlockSize
          daux=COS(Stack(iBlock,StackPtr))
          IF (daux == 0._DP) THEN
            EvalErrType=1
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr)=1._DP/daux
          END IF
        END DO
!$omp end parallel do

        
      CASE (cSin)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=SIN(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do

        
      CASE(cSinh)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=SINH(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
      

      CASE(cSqrt)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
          IF (Stack(iBlock,StackPtr) < 0._DP) THEN
            EvalErrType=3
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr)=SQRT(Stack(iBlock,StackPtr))
          END IF
        END DO
!$omp end parallel do


      CASE (cTan)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=TAN(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do

        
      CASE (cTanh)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=TANH(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        

        !------------------------------------------------------------
        ! Misc
        !------------------------------------------------------------
      CASE (cImmed)
        StackPtr=StackPtr+1
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=Comp%Immed(DataPtr)
        END DO
!$omp end parallel do
        DataPtr=DataPtr+1


      CASE (cJump)
        DataPtr=Comp%ByteCode(InstPtr+2)
        InstPtr=Comp%ByteCode(InstPtr+1)


        !------------------------------------------------------------
        ! Operators
        !------------------------------------------------------------
      CASE (cNeg)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=-Stack(iBlock,StackPtr)
        END DO
!$omp end parallel do


      CASE (cAdd)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)+Stack(iBlock,StackPtr)
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1
        

      CASE (cSub)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)-Stack(iBlock,StackPtr)
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1


      CASE (cMul)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)*Stack(iBlock,StackPtr)
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1
        

      CASE (cDiv)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
          IF (Stack(iBlock,StackPtr) == 0._DP) THEN
            EvalErrType=1
            Res(iBlock)=zero
          ELSE
            Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)/Stack(iBlock,StackPtr)
          END IF
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1
        

      CASE (cMod)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
          IF (Stack(iBlock,StackPtr) == 0._DP) THEN
            EvalErrType=1
            Res(iBlock)=Zero
          ELSE
            Stack(iBlock,StackPtr-1)=MOD(Stack(iBlock,StackPtr-1),Stack(iBlock,StackPtr))
          END IF
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1


      CASE (cPow)
!$omp parallel do default(shared) private(iBlock)
        DO  iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=Stack(iBlock,StackPtr-1)**Stack(iBlock,StackPtr)
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1

        
      CASE (cEqual)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) == Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1


      CASE (cNEqual)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) /= Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1


      CASE (cLess)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) < Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1


      CASE (cLessOrEq)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) <= Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1

        
      CASE (cGreater)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) > Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1

        
      CASE (cGreaterOrEq)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(Stack(iBlock,StackPtr-1) >= Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1

        
      CASE (cAnd)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(DbleToLogc(Stack(iBlock,StackPtr-1)) .AND. &
                DbleToLogc(Stack(iBlock,StackPtr)) )
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1


      CASE (cOr)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr-1)=LogcToDble(DbleToLogc(Stack(iBlock,StackPtr-1)) .OR. &
                DbleToLogc(Stack(iBlock,StackPtr)) )
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1


      CASE (cNot)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=LogcToDble( .NOT. DbleToLogc(Stack(iBlock,StackPtr)) )
        END DO
!$omp end parallel do

        
        !------------------------------------------------------------
        ! Degrees-radians conversion
        !------------------------------------------------------------
      CASE (cDeg)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=RadToDeg(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do

        
      CASE (cRad)
!$omp parallel do default(shared) private(iBlock)
        DO iBlock=1,BlockSize
           Stack(iBlock,StackPtr)=DegToRad(Stack(iBlock,StackPtr))
        END DO
!$omp end parallel do
        

      CASE DEFAULT
        StackPtr  = StackPtr+1
        iVariable = Comp%ByteCode(InstPtr)-VarBegin+1

        ! Do we have to process one of the scalar variables of one of the block variables
        IF (iVariable >= istartValScalar) THEN
!$omp parallel do default(shared) private(iBlock)
          DO iBlock=1,BlockSize
             Stack(iBlock,StackPtr)=ValScalar(iVariable-istartValScalar+1)
          END DO
!$omp end parallel do
        ELSE
          IF (iDim == 1) THEN
!$omp parallel do default(shared) private(iBlock)
            DO iBlock=1,BlockSize
               Stack(iBlock,StackPtr)=ValBlock(iBlock,iVariable)
            END DO
!$omp end parallel do
          ELSE
!$omp parallel do default(shared) private(iBlock)
            DO iBlock=1,BlockSize
               Stack(iBlock,StackPtr)=ValBlock(iVariable,iBlock)
            END DO
!$omp end parallel do
          END IF
        END IF
      END SELECT
    END DO
    
    EvalErrType = 0
!$omp parallel do default(shared) private(iBlock)
    DO iBlock=1,BlockSize
       Res(iBlock) = Stack(iBlock,StackPtr)
    END DO
!$omp end parallel do
  END SUBROUTINE evalFunctionBlock
  
  ! *****************************************************************************

!<function>

  FUNCTION fparser_ErrorMsg (EvalErrType) RESULT (msg)

!<description>
    ! Return error message of function parser
!</description>

    ! local constants
    CHARACTER (LEN=*), DIMENSION(6), PARAMETER :: m = (/ 'Division by zero                  ', &
                                                         'Argument of SQRT negative         ', &
                                                         'Argument of LOG negative          ', &
                                                         'Argument of ASIN or ACOS illegal  ', &
                                                         'Argument of ASINH or ACOSH illegal', &
                                                         'Argument of ATANH illegal         ' /)

!<input>
    ! Error identifier
    INTEGER, INTENT(IN) :: EvalErrType
!</input>   
    
!<result>
    ! Error messages
    CHARACTER (LEN=LEN(m)) :: msg
!</result>
!</function>

    IF (EvalErrType < 1 .OR. EvalErrType > SIZE(m)) THEN
      msg = ''
    ELSE
      msg = m(EvalErrType)
    ENDIF
  END FUNCTION fparser_ErrorMsg

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_PrintByteCode(rparser, iComp)

!<description>
    ! Print the compiled bytecode stack
!</description>

!<input>
    ! Function parser
    TYPE(t_fparser), INTENT(IN) :: rparser

    ! Function identifier
    INTEGER, INTENT(IN) :: iComp
!</input>
!</subroutine>

    ! local variables
    TYPE(t_fparserComponent), POINTER :: Comp
    CHARACTER(LEN=5) :: n
    INTEGER :: InstPtr, DataPtr, StackPtr,nparams
    INTEGER(is) :: opCode
    
    nparams  = 1
    DataPtr  = 1
    StackPtr = 0
    InstPtr  = 0
    Comp     => rparser%Comp(iComp)

    DO WHILE(InstPtr < Comp%ByteCodeSize)
      InstPtr=InstPtr+1
      
      WRITE(*,FMT='(I8.8,1X,":",1X)',ADVANCE="NO") InstPtr
      
      opCode=Comp%ByteCode(InstPtr)
      SELECT CASE(opCode)
        
      CASE(cIf)
        WRITE(*,FMT='(A,1X,T10,I8.8)') "jz",Comp%ByteCode(InstPtr+1)+1
        InstPtr=InstPtr+2

      CASE(cJump)
        WRITE(*,FMT='(A,1X,T10,I8.8)') "jump",Comp%ByteCode(InstPtr+1)+1
        InstPtr=InstPtr+2

      CASE(cImmed)
        WRITE(*,FMT='(A,1X,T10,G8.2)') "push",Comp%Immed(DataPtr)
        DataPtr=DataPtr+1

      CASE DEFAULT
        IF (opCode < VarBegin) THEN
          SELECT CASE(opCode)
          CASE(cNEG); n="neg"
          CASE(cADD); n="add"
          CASE(cSUB); n="sub"
          CASE(cMUL); n="mul"
          CASE(cDIV); n="div"
          CASE(cMOD); n="mod"
          CASE(cPOW); n="pow"
          CASE(cEqual); n="eq"
          CASE(cNEqual); n="ne"
          CASE(cLess); n="lt"
          CASE(cLessOrEq); n="le"
          CASE(cGreater); n="gt"
          CASE(cGreaterOrEq); n="ge"
          CASE(cAND); n="and"
          CASE (cOR); n="or"
          CASE(cNOT); n="not"
          CASE(cDEG); n="deg"
          CASE(cRAD); n="rad"

          CASE DEFAULT
            n       = Funcs(opCode)
            nparams = MathFunctionParameters(opCode)
          END SELECT
          WRITE(*,FMT='(A,T10,A,"  (",I1,") ")') TRIM(n),"Par",nparams
          
        ELSE
          WRITE(*,FMT='(A,T10,A,1X,I4.4)') "push","Var",opCode-VarBegin+1
        END IF
        
      END SELECT
    END DO
  END SUBROUTINE fparser_PrintByteCode

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE CheckSyntax (Func,Var)

!<description>
    ! Check syntax of function string, returns 0 if syntax is ok
!</description>
    
!<input>
    ! Function string without spaces
    CHARACTER (LEN=*), INTENT(in) :: Func

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>
!</subroutine>
    
    ! local avriables
    TYPE(t_stack) :: functionParenthDepth
    INTEGER(is) :: n,opSize
    CHARACTER (LEN=1) :: c
    REAL(DP) :: r
    LOGICAL :: err
    INTEGER :: FuncIndex,FuncIndex2
    INTEGER :: ParCnt,ib,in,FuncLength,idummy
    
    ! Initialization
    FuncIndex = 1
    ParCnt = 0
    FuncLength = LEN_TRIM(Func)
    CALL stack_create(functionParenthDepth,MAX(5,INT(FuncLength/4._DP)),ST_INT)

    DO
      IF (FuncIndex > FuncLength) THEN
        PRINT *, "CheckSyntax: invalid function string!"
        CALL sys_halt()
      END IF
      c = Func(FuncIndex:FuncIndex)

      ! Check for valid operand (must appear)

      ! Check for leading - or !
      IF (c == '-' .OR. c == '!') THEN                      
        FuncIndex = FuncIndex+1
        IF (FuncIndex > FuncLength) THEN
          PRINT *, "CheckSyntax: Premature end of string!"
          CALL sys_halt()
        END IF
        c = Func(FuncIndex:FuncIndex)
      END IF
            
      ! Check for math function
      n = MathFunctionIndex (Func(FuncIndex:))
      IF (n > 0) THEN
        ! Math function found
        FuncIndex = FuncIndex+LEN_TRIM(Funcs(n))
        IF (FuncIndex > FuncLength) THEN
          PRINT *, "CheckSyntax: Premature end of string!"
          CALL sys_halt()
        END IF
        c = Func(FuncIndex:FuncIndex)
        IF (c /= '(') THEN
          PRINT *, "CheckSyntax: Expecting ( after function!",Func(FuncIndex:)
          CALL sys_halt()
        END IF
        FuncIndex2=FuncIndex+1
        IF (FuncIndex2 > FuncLength) THEN
          PRINT *, "CheckSyntax: Premature end of string!"
          CALL sys_halt()
        END IF
        IF (Func(FuncIndex2:FuncIndex2) == ')') THEN
          FuncIndex=FuncIndex2+1
          IF (FuncIndex > FuncLength) THEN
            PRINT *, "CheckSyntax: Premature end of string!"
            CALL sys_halt()
          END IF
          c = Func(FuncIndex:FuncIndex)
          
          ! Ugly, but other methods would just be uglier ...
          GOTO 999
        END IF
        
        ! Push counter for parenthesss to stack
        CALL stack_pushback(functionParenthDepth,ParCnt+1)
      END IF
      
      ! Check for opening parenthesis
      IF (c == '(') THEN
        ParCnt = ParCnt+1
        FuncIndex = FuncIndex+1
        IF (FuncIndex > FuncLength) THEN
          PRINT *, "CheckSyntax: Premature end of string!"
          CALL sys_halt()
        END IF
        IF (Func(FuncIndex:FuncIndex) == ')') THEN
          PRINT *, "CheckSyntax: Empty parantheses!"
          CALL sys_halt()
        END IF
        CYCLE
      END IF

      ! Check for number
      IF (SCAN(c,'0123456789.') > 0) THEN
        r = RealNum (Func(FuncIndex:),ib,in,err)
        IF (err) THEN
          PRINT *, "CheckSyntax: Invalid number format!"
          CALL sys_halt()
        END IF
        FuncIndex = FuncIndex+in-1
        IF (FuncIndex > FuncLength) EXIT
        c = Func(FuncIndex:FuncIndex)
      ELSEIF (c == '_') THEN
        ! Check for constant
        n = ConstantIndex (Func(FuncIndex:))
        IF (n == 0) THEN
          PRINT *, "CheckSyntax: Invalid constant!"
          CALL sys_halt()
        END IF
        FuncIndex = FuncIndex+LEN_TRIM(Consts(n))+1
        IF (FuncIndex > FuncLength) EXIT
        c=Func(FuncIndex:FuncIndex)
      ELSEIF (c == '@') THEN
        ! Check for expression
        n = ExpressionIndex (Func(FuncIndex:))
        IF (n == 0) THEN
          PRINT *, "CheckSyntax: Invalid expression!"
          CALL sys_halt()
        END IF
        CALL CheckSyntax(ExpressionStrings(n), Var)
        FuncIndex = FuncIndex+LEN_TRIM(Expressions(n))+1
        IF (FuncIndex > FuncLength) EXIT
        c=Func(FuncIndex:FuncIndex)
      ELSE
        ! Check for variable
        n = VariableIndex (Func(FuncIndex:),Var,ib,in)
        IF (n == 0) THEN
          PRINT *, "CheckSyntax: Invalid element!"
          CALL sys_halt()
        END IF
        FuncIndex = FuncIndex+in-1
        IF (FuncIndex > FuncLength) EXIT
        c = Func(FuncIndex:FuncIndex)
      END IF

      ! Check for closing parenthesis
      DO WHILE (c == ')')
        IF (.NOT.stack_isempty(functionParenthDepth)) THEN
          IF(stack_backInt(functionParenthDepth) == ParCnt) &
              idummy=stack_popbackInt(functionParenthDepth)
        END IF
        ParCnt = ParCnt-1
        IF (ParCnt < 0) THEN
          PRINT *, "CheckSyntax: Mismatched parenthesis!"
          CALL sys_halt()
        END IF
        IF (Func(FuncIndex-1:FuncIndex-1) == '(') THEN
          PRINT *, "CheckSyntax: Empty parentheses!"
          CALL sys_halt()
        END IF
        FuncIndex = FuncIndex+1
        IF (FuncIndex > FuncLength) EXIT
        c = Func(FuncIndex:FuncIndex)
      END DO

      ! Now, we have a legal operand: A legal operator or end of string must follow
999   IF (FuncIndex > FuncLength) EXIT
      
      ! Check operators
      opSize = 0
      IF (.NOT.stack_isempty(functionParenthDepth)) THEN
        IF (c == ',' .AND. &
            stack_backInt(functionParenthDepth) == parCnt) THEN
          opSize = 1
        ELSE
          opSize = isOperator(Func(FuncIndex:))
        END IF
      ELSE
        opSize = isOperator(Func(FuncIndex:))
      END IF
      IF (opSize == 0) THEN
        PRINT *, "CheckSyntax: Operator expected!"
        CALL sys_halt()
      END IF
      
      ! Now, we have an operand and an operator: the next loop will check for another 
      ! operand (must appear)
      FuncIndex = FuncIndex+opSize
    END DO
    IF (ParCnt > 0) THEN
      PRINT *, "CheckSyntax: Missing )!"
      CALL sys_halt()
    END IF
    CALL stack_release(functionParenthDepth)
  END SUBROUTINE CheckSyntax

  ! *****************************************************************************

!<function>

  FUNCTION isOperator (str) RESULT (n)

!<description>
    ! Return 0 if given string is not an operator, else the size of the
    ! operator
!</description>
    
!<input>
    ! Operator string
    CHARACTER(LEN=*), INTENT(IN) :: str    
!</input>

!<result>
    ! Length of operator, 0 if string is no operator
    INTEGER(is) :: n
!</result>
!</function>

    ! local variables
    INTEGER(is) :: j,m

    n = 0
    DO j=cAdd,cOr
      m = LEN_TRIM(Ops(j))
      IF (str(1:m) == TRIM(Ops(j))) THEN
        n = m
        EXIT
      END IF
    END DO
    
  END FUNCTION isOperator

  ! *****************************************************************************

!<function>

  FUNCTION MathFunctionIndex (str) RESULT (n)

!<description>
    ! Return index of math function beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    CHARACTER (LEN=*), INTENT(in) :: str
!</input>

!<result>
    ! Index of math function
    INTEGER(is) :: n
!</result>
!</function>

    ! local variables
    INTEGER(is) :: j
    INTEGER k
    CHARACTER (LEN=LEN(Funcs)) :: fun
    
    ! Check all math functions
    n = 0
    DO j=cIf,cSec
      k = MIN(LEN_TRIM(Funcs(j)), LEN(str))

      ! Compare lower case letters
      CALL sys_tolower(str(1:k), fun)
      IF (fun == Funcs(j)) THEN
        ! Found a matching function
        n = j
        RETURN
      END IF
    END DO
  END FUNCTION MathFunctionIndex

  ! *****************************************************************************
  
!<function>

  FUNCTION MathFunctionParameters (FuncIndex) RESULT (nparameters)

!<description>
    ! Return number of required parameters
!</description>

!<input>
    ! Index of function
    INTEGER(is) :: FuncIndex
!</input>

!<result>
    ! Number if required parameters
    INTEGER :: nparameters
!</result>
!</function>

    SELECT CASE(FuncIndex)
    CASE(cIf)
      nparameters = 3
    CASE(cMin:cAtan2)
      nparameters = 2
    CASE(cAbs:cSec)
      nparameters = 1
    CASE DEFAULT
      nparameters=0
      PRINT *, "*** MathFunctionParamters: Not a function"
    END SELECT
  END FUNCTION MathFunctionParameters

  ! *****************************************************************************

!<function>

  FUNCTION ConstantIndex (str) RESULT (n)

!<description>
    ! Return index of predefined constants beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    CHARACTER (LEN=*), INTENT(in) :: str
!</input>

!<result>
    ! Index of math function
    INTEGER(is) :: n
!</result>
!</function>

    ! local variables
    INTEGER(is) :: j
    INTEGER k
    CHARACTER (LEN=LEN(Consts)) :: con
    
    ! Check all predefined constants
    n = 0
    DO j=1,nConsts
      k = MIN(LEN_TRIM(Consts(j)), LEN(str(2:)))

      ! Compare lower case letters
      CALL sys_tolower(str(2:k+1),con)
      IF (con == Consts(j)) THEN
        ! Found a matching constant
        n = j
        RETURN
      END IF
    END DO
  END FUNCTION ConstantIndex

  ! *****************************************************************************

!<function>

  FUNCTION ExpressionIndex (str) RESULT (n)

!<description>
    ! Return index of predefined expression beginnig at 1st position of string str
!</description>

!<input>
    ! Math function string
    CHARACTER (LEN=*), INTENT(in) :: str
!</input>

!<result>
    ! Index of math function
    INTEGER(is) :: n
!</result>
!</function>

    ! local variables
    INTEGER(is) :: j
    INTEGER k
    CHARACTER (LEN=LEN(Expressions)) :: expression

    ! Check all predefined expressions
    n = 0
    DO j=1,nExpressions
      k = MIN(LEN_TRIM(Expressions(j)), LEN(str(2:)))

      ! Compare lower case letters
      CALL sys_tolower(str(2:k+1),expression)

      IF (expression == Expressions(j)) THEN
        ! Found a matching expression
        n = j
        RETURN
      END IF
    END DO
  END FUNCTION ExpressionIndex

  ! *****************************************************************************

!<function>
  
  FUNCTION VariableIndex (str, Var, ibegin, inext) RESULT (n)

!<description>
    ! Return index of variable at begin of string str (returns 0 if no variable found)
!</description>

!<input>
    ! String
    CHARACTER (LEN=*), INTENT(in) :: str

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>

!<output>
    ! OPTIONAL: Start position of variable name
    INTEGER, OPTIONAL, INTENT(out) :: ibegin

    ! OPTIONAL: Position of character after name
    INTEGER, OPTIONAL, INTENT(out) :: inext
!</output>

!<result>
    ! Index of variable
    INTEGER(is) :: n
!</result>
!</function>

    ! local variables    
    INTEGER :: j,ib,in,lstr
    
    n = 0
    lstr = LEN_TRIM(str)
    IF (lstr > 0) THEN
      ! Search for first character in str
      DO ib=1,lstr
        ! When lstr>0 at least 1 char in str
        IF (str(ib:ib) /= ' ') EXIT
      END DO
      
      ! Search for name terminators
      DO in=ib,lstr
        IF (SCAN(str(in:in),'+-*/%^),&|<>=! ') > 0) EXIT
      END DO
      DO j=1,SIZE(Var)
        IF (str(ib:in-1) == Var(j)) THEN
          ! Variable name found
          n = j
          EXIT
        END IF
      END DO
    END IF

    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
  END FUNCTION VariableIndex

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE RemoveSpaces (str)

!<description>
    ! Remove Spaces from string, remember positions of characters in
    ! old string
!</description>
   
!<inputoutput>
    ! String from which spaces should be removed
    CHARACTER (LEN=*), INTENT(inout) :: str
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: k,lstr
    
    lstr = LEN_TRIM(str)
    k = 1
    DO WHILE (str(k:lstr) /= ' ')                             
      IF (str(k:k) == ' ') THEN
        str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
        k = k-1
      END IF
      k = k+1
    END DO
  END SUBROUTINE RemoveSpaces

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE Replace (ca,cb,str)

!<description>
    ! Replace ALL appearances of character set ca in string str by character set cb
!</description>

!<input>
    ! Source characters
    CHARACTER (LEN=*), INTENT(in) :: ca
    
    ! Destination characters
    CHARACTER (LEN=LEN(ca)), INTENT(in) :: cb
!</input>

!<inputoutput>
    ! String
    CHARACTER (LEN=*), INTENT(inout) :: str
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER :: j,lca
    
    lca = LEN(ca)
    DO j=1,LEN_TRIM(str)-lca+1
      IF (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    END DO
  END SUBROUTINE Replace

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE Compile (Comp, Func, Var)

!<description>
    ! Compile i-th function string Func into bytecode
!</description>

!<input>
    ! Function string
    CHARACTER (LEN=*), INTENT(in) :: Func

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>

!<inputoutput>
    ! Function parser component
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(is), DIMENSION(:), POINTER :: ByteCode
    REAL(DP), DIMENSION(:), POINTER :: Immed
    INTEGER :: ind,isize
    
    ! (Re-)initialize the bytecode structure (if required)
    IF (ASSOCIATED(Comp%ByteCode)) DEALLOCATE (Comp%ByteCode)
    IF (ASSOCIATED(Comp%Immed))    DEALLOCATE (Comp%Immed)
    Comp%ByteCodeSize = 0
    Comp%ImmedSize    = 0
    Comp%StackSize    = 0
    Comp%StackPtr     = 0
    Comp%isVectorizable = .TRUE.

    ! Neither the stack for the bytecode nor the stack for the immediate expressions
    ! can exceed the size of the function string. Hence, allocate some initial memory
    isize=FunctionSize(Func)
    ALLOCATE(Comp%ByteCode(isize),Comp%Immed(isize))

    ! Compile function string into bytecode
    ind = CompileExpression(Comp,Func,1,Var)
    
    ! Adjust memory size of bytecode stack
    ALLOCATE(ByteCode(Comp%ByteCodeSize))
    ByteCode=Comp%ByteCode
    DEALLOCATE(Comp%ByteCode)
    ALLOCATE(Comp%ByteCode(Comp%ByteCodeSize))
    Comp%ByteCode=ByteCode
    DEALLOCATE(ByteCode)
    
    ! Adjust memory size of immediate stack
    ALLOCATE(Immed(Comp%ImmedSize))
    Immed=Comp%Immed
    DEALLOCATE(Comp%Immed)
    ALLOCATE(Comp%Immed(Comp%ImmedSize))
    Comp%Immed=Immed
    DEALLOCATE(Immed)
  END SUBROUTINE Compile

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE incStackPtr (Comp)

!<description>
    ! Increase stack pointer
!</description>

!<inputoutput>
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    Comp%StackPtr=Comp%StackPtr+1
    IF (Comp%StackPtr > Comp%StackSize) Comp%StackSize=Comp%StackSize+1
  END SUBROUTINE incStackPtr

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE AddCompiledByte (Comp, byte)

!<description>
    ! Add compiled byte to bytecode
!</description>

!<input>
    ! Value of byte to be added
    INTEGER(is), INTENT(in) :: byte
!</input>

!<inputoutput>
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP) :: daux

    Comp%ByteCodeSize = Comp%ByteCodeSize + 1    
    Comp%ByteCode(Comp%ByteCodeSize) = byte

    ! Try to optimize the compiled bytecode. Check the bytecode instruction and
    ! compute some values on-the-fly of this is possible
    SELECT CASE(byte)
      !------------------------------------------------------------
      ! Functions
      !------------------------------------------------------------
    CASE (cAbs)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=ABS(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cAcos)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) < -1._DP .OR. &
            Comp%Immed(Comp%ImmedSize) >  1._DP) THEN
          PRINT *, "*** Parser error: invalid argument for ACOS!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=ACOS(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cAsin)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) < -1._DP .OR. &
            Comp%Immed(Comp%ImmedSize) >  1._DP) THEN
          PRINT *, "*** Parser error: invalid argument for ACOS!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=ASIN(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cAtan)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=ATAN(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cAtan2)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=ATAN2(Comp%Immed(Comp%ImmedSize),Comp%Immed(Comp%ImmedSize-1))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cAcosh)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        daux=Comp%Immed(Comp%ImmedSize)+SQRT(Comp%Immed(Comp%ImmedSize)**2-1)
        IF (daux <= 0) THEN
          PRINT *, "*** Parser error: invalid argument for ACOSH!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=LOG(daux)
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cAnint)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=ANINT(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cAint)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=AINT(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cAsinh)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        daux=Comp%Immed(Comp%ImmedSize)+SQRT(Comp%Immed(Comp%ImmedSize)**2-1)
        IF (daux <= 0) THEN
          PRINT *, "*** Parser error: invalid argument for ASINH!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=LOG(daux)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cAtanh)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) == -1._DP) THEN
          PRINT *, "*** Parser error: invalid argument for ATANH!"
          CALL sys_halt()
        END IF
        daux=(1+Comp%Immed(Comp%ImmedSize))/(1-Comp%Immed(Comp%ImmedSize))
        IF (daux <= 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for ATANH!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=LOG(daux)/2._DP
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cCeil)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=CEILING(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cCos)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=COS(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cCosh)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=COSH(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cCot)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        daux=TAN(Comp%Immed(Comp%ImmedSize))
        IF (daux == 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for COT!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=1/daux
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cCsc)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        daux=SIN(Comp%Immed(Comp%ImmedSize))
        IF (daux == 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for CSC!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=1/daux
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cExp)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=EXP(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cFloor)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=FLOOR(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cIf)
      ! No optimization possible

    CASE (cLog)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) <= 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for LOG!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=LOG(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cLog10)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) <= 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for LOG!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=LOG10(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cMax)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=MAX(Comp%Immed(Comp%ImmedSize),Comp%Immed(Comp%ImmedSize-1))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cMin)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=MIN(Comp%Immed(Comp%ImmedSize),Comp%Immed(Comp%ImmedSize-1))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cSec)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        daux=COS(Comp%Immed(Comp%ImmedSize))
        IF (daux == 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for SEC!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=1/daux
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cSin)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=SIN(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cSinh)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=SINH(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cSqrt)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) < 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for SQRT!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize)=SQRT(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cTan)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=TAN(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cTanh)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=TANH(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

      !------------------------------------------------------------
      ! Misc
      !------------------------------------------------------------
    CASE (cImmed,cJump)
      ! No optimization needed

      !------------------------------------------------------------
      ! Operators
      !------------------------------------------------------------
    CASE (cNeg)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=-(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cAdd)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)+Comp%Immed(Comp%ImmedSize)
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cSub)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)-Comp%Immed(Comp%ImmedSize)
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cMul)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)*Comp%Immed(Comp%ImmedSize)
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cDiv)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) == 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for DIV!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)/Comp%Immed(Comp%ImmedSize)
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cMod)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        IF (Comp%Immed(Comp%ImmedSize) == 0._DP) THEN
          PRINT *, "*** Parser error: invalid argument for MOD!"
          CALL sys_halt()
        END IF
        Comp%Immed(Comp%ImmedSize-1)=MOD(Comp%Immed(Comp%ImmedSize-1),Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cPow)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=Comp%Immed(Comp%ImmedSize-1)**Comp%Immed(Comp%ImmedSize)
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cEqual)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) == Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cNEqual)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) /= Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cLess)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) < Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cLessOrEq)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) <= Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cGreater)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) > Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cGreaterOrEq)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(Comp%Immed(Comp%ImmedSize-1) >= Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF
      
    CASE (cAnd)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(DbleToLogc(Comp%Immed(Comp%ImmedSize-1)) .AND.&
            DbleToLogc(Comp%Immed(Comp%ImmedSize)))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cOr)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed .AND.&
          Comp%ByteCode(Comp%ByteCodeSize-2) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize-1)=LogcToDble(DbleToLogc(Comp%Immed(Comp%ImmedSize-1)) .OR.&
            DbleToLogc(Comp%Immed(Comp%ImmedSize)))
        CALL RemoveCompiledImmediate(Comp)
        CALL RemoveCompiledByte(Comp)
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cNot)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=LogcToDble(.NOT.DbleToLogc(Comp%Immed(Comp%ImmedSize)))
        CALL RemoveCompiledByte(Comp)
      END IF

      !------------------------------------------------------------
      ! Degrees-radians conversion
      !------------------------------------------------------------
    CASE (cDeg)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=RadToDeg(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF

    CASE (cRad)
      IF (Comp%ByteCode(Comp%ByteCodeSize-1) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize)=DegToRad(Comp%Immed(Comp%ImmedSize))
        CALL RemoveCompiledByte(Comp)
      END IF      
    END SELECT
  END SUBROUTINE AddCompiledByte

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE RemoveCompiledByte (Comp)

!<description>
    ! Remove last compiled byte from bytecode
!</description>

!<inputoutput>
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>
   
    Comp%ByteCode(Comp%ByteCodeSize) = 0
    Comp%ByteCodeSize = Comp%ByteCodeSize - 1
  END SUBROUTINE RemoveCompiledByte

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE AddImmediate (Comp, immediate)

!<description>
    ! Add immediate
!</description>

!<input>
    ! Value of byte to be added
    REAL(DP), INTENT(in) :: immediate
!</input>

!<inputoutput>
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>
       
    Comp%ImmedSize = Comp%ImmedSize + 1
    Comp%Immed(Comp%ImmedSize) = immediate
  END SUBROUTINE AddImmediate

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE RemoveCompiledImmediate (Comp)

!<description>
    ! Remove last compiled immediate from immediate stack
!</description>

!<inputoutput>
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    Comp%Immed(Comp%ImmedSize) = 0
    Comp%ImmedSize = Comp%ImmedSize - 1
  END SUBROUTINE RemoveCompiledImmediate

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE AddFunctionOpcode (Comp, opcode)

!<description>
    ! Add compiled byte to bytecode
!</description>

!<input>
    ! Value of opcode to be added
    INTEGER(is), INTENT(in) :: opcode
!</input>

!<inputoutput>
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>
    
    IF (Comp%useDegreeConversion) THEN
      SELECT CASE(opCode)
      CASE(cCos,Ccosh,cCot,cCsc,cSec,cSin,cSinh,cTan,cTanh)
        CALL AddCompiledByte(Comp,cRad)
      END SELECT
    END IF
    
    CALL AddCompiledByte(Comp,opcode)
    
    IF (Comp%useDegreeConversion) THEN
      SELECT CASE(opCode)
      CASE(cAcos,cAcosh,cAsinh,cAtanh,cAsin,cAtan,cAtan2)
        CALL AddCompiledByte(Comp, cDeg)
      END SELECT
    END IF
  END SUBROUTINE AddFunctionOpcode

  ! *****************************************************************************

!<function>

  FUNCTION RealNum (str, ibegin, inext, error) RESULT (res)

!<description>
    ! Get real number from string
    ! Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
!</description>

!<input>
    ! String
    CHARACTER (LEN=*), INTENT(in) :: str
!</input>

!<output>
    ! OPTIONAL: Start position of real number
    INTEGER, OPTIONAL, INTENT(out) :: ibegin
    
    ! OPTIONAL: 1st character after real number
    INTEGER, OPTIONAL, INTENT(out) :: inext

    ! OPTIONAL: Error flag
    LOGICAL, OPTIONAL, INTENT(out) :: error
!</output>

!<result>
    ! Real number
    REAL(DP) :: res
!</result>
!</function>
 
    ! local variables
    INTEGER :: ib,in,istat
    LOGICAL :: Bflag,               & ! .T. at begin of number in str
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
    DO WHILE (in <= LEN_TRIM(str))
      SELECT CASE (str(in:in))
      CASE (' ') ! Only leading blanks permitted
        ib = ib+1
        IF (InMan .OR. Eflag .OR. InExp) EXIT
      CASE ('+','-') ! Permitted only
        IF     (Bflag) THEN           
          InMan=.true.; Bflag=.false. ! - at beginning of mantissa
        ELSEIF (Eflag) THEN               
          InExp=.true.; Eflag=.false. ! - at beginning of exponent
        ELSE
          EXIT ! - otherwise CALL sys_halt()
        ENDIF
      CASE ('0':'9') ! Mark
        IF     (Bflag) THEN           
          InMan=.true.; Bflag=.false. ! - beginning of mantissa
        ELSEIF (Eflag) THEN               
          InExp=.true.; Eflag=.false. ! - beginning of exponent
        ENDIF
        IF (InMan) DInMan=.true. ! Mantissa contains digit
        IF (InExp) DInExp=.true. ! Exponent contains digit
      CASE ('.')
        IF     (Bflag) THEN
          Pflag=.true. ! - mark 1st appearance of '.'
          InMan=.true.; Bflag=.false. !   mark beginning of mantissa
        ELSEIF (InMan .AND..NOT.Pflag) THEN
          Pflag=.true. ! - mark 1st appearance of '.'
        ELSE
          EXIT ! - otherwise CALL sys_halt()
        END IF
      CASE ('e','E','d','D') ! Permitted only
        IF (InMan) THEN
          Eflag=.true.; InMan=.false. ! - following mantissa
        ELSE
          EXIT ! - otherwise CALL sys_halt()
        ENDIF
      CASE DEFAULT
        EXIT ! CALL sys_halt() at all other characters
      END SELECT
      in = in+1
    END DO
    err = (ib > in-1) .OR. (.NOT.DInMan) .OR.&
        &((Eflag.OR.InExp).AND..NOT.DInExp)
    IF (err) THEN
      res = 0.0_DP
    ELSE
      READ(str(ib:in-1),*,IOSTAT=istat) res
      err = istat /= 0
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
    IF (PRESENT(error))  error  = err
  END FUNCTION RealNum

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION FunctionSize (Func) RESULT (fsize)

!<description>
    ! Return the size of the total function including external expressions
!</description>
    
!<input>
    ! Function string
    CHARACTER (LEN=*),  INTENT(in) :: Func
!</input>

!<result>
    ! Size of function string
    INTEGER :: fsize
!</result>
!</function>
    
    ! local variables
    INTEGER :: ind,n
    CHARACTER(LEN=1) :: c

    ! Determine size of given expression
    fsize=LEN_TRIM(Func)

    ! "Parse" string for externally defined expressions
    DO ind=1,fsize
      c = Func(ind:ind)
      IF (c == '@') THEN
        n = ExpressionIndex (Func(ind:))
        IF (n == 0) THEN
          PRINT *, "FunctionSize: Invalid expression!"
          CALL sys_halt()
        END IF
        fsize = fsize+FunctionSize(ExpressionStrings(n))
      END IF
    END DO
  END FUNCTION FunctionSize

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileExpression(Comp, Func, ind, Var, stopAtComma) RESULT(ind2)

!<description>
    ! Compiles ','
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var

    ! OPTIONAL: stop at comma
    LOGICAL, INTENT(IN), OPTIONAL :: stopAtComma
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    INTEGER :: FuncLength
    
    FuncLength=LEN_TRIM(Func)
    
    ind2 = CompileOr(Comp,Func, ind, Var)
    IF(ind2 > FuncLength) RETURN

    IF (PRESENT(stopAtComma)) THEN
      IF (stopAtComma) RETURN
    END IF
    
    DO WHILE (Func(ind2:ind2) == ',')
      ind2 = CompileOr(Comp, Func, ind2+1, Var)
      IF (ind2 > FuncLength) RETURN
    END DO
  END FUNCTION CompileExpression

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileOr(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles '|'
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    INTEGER :: FuncLength

    FuncLength=LEN_TRIM(Func)

    ind2 = CompileAnd(Comp, Func, ind, Var)
    IF (ind2 > FuncLength) RETURN

    DO WHILE(Func(ind2:ind2) == '|')
      ind2 = CompileAnd(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, cOr)
      Comp%StackPtr = Comp%StackPtr-1
      IF (ind2 > FuncLength) RETURN
    END DO
  END FUNCTION CompileOr

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileAnd(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles '&'
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    INTEGER :: FuncLength

    FuncLength=LEN_TRIM(Func)
    
    ind2 = CompileComparison(Comp, Func, ind, Var)
    IF (ind2 > FuncLength) RETURN
    
    DO WHILE(Func(ind2:ind2) == '&')
      ind2 = CompileComparison(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, cAnd)
      Comp%StackPtr = Comp%StackPtr-1
      IF (ind2 > FuncLength) RETURN
    END DO
  END FUNCTION CompileAnd

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileComparison(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles '=', '<' and '>'
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    CHARACTER(LEN=1) :: c
    INTEGER(is) :: opSize
    INTEGER :: FuncLength

    FuncLength=LEN_TRIM(Func)
    
    ind2 = CompileAddition(Comp, Func, ind, Var)
    IF (ind2 > FuncLength) RETURN
    
    c=Func(ind2:ind2)
    DO WHILE(c == '=' .OR. c == '<' .OR. c == '>' .OR. c == '!')
      opSize = MERGE(2,1,Func(ind2+1:ind2+1) == '=')
      ind2 = CompileAddition(Comp, Func, ind2+opSize, Var)

      SELECT CASE(c)
      CASE('=')
        CALL AddCompiledByte(Comp, cEqual)

      CASE('<')
        CALL AddCompiledByte(Comp, MERGE(cLess,cLessOrEq,opSize==1))
        
      CASE('>')
        CALL AddCompiledByte(Comp, MERGE(cGreater,cGreaterOrEq,opSize==1))
        
      CASE('!')
        CALL AddCompiledByte(Comp, cNEqual)
      END SELECT
      Comp%StackPtr = Comp%StackPtr-1
      
      IF (ind2 > FuncLength) RETURN
      c=Func(ind2:ind2)
    END DO
  END FUNCTION CompileComparison

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileAddition(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles '+' and '-'
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    CHARACTER(LEN=1) :: c
    INTEGER :: FuncLength

    FuncLength=LEN_TRIM(Func)
    
    ind2 = CompileMult(Comp, Func, ind, Var)
    IF (ind2 > FuncLength) RETURN
    
    c=Func(ind2:ind2)
    DO WHILE(c == '+' .OR. c == '-')
      ind2 = CompileMult(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, MERGE(cAdd,cSub,c == '+'))
      Comp%StackPtr = Comp%StackPtr-1
      
      IF (ind2 > FuncLength) RETURN
      c=Func(ind2:ind2)
    END DO
  END FUNCTION CompileAddition

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileMult(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles '*', '/' and '%'
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    CHARACTER(LEN=1) :: c
    INTEGER :: FuncLength

    FuncLength=LEN_TRIM(Func)
    
    ind2 = CompileUnaryMinus(Comp, Func, ind, Var)
    IF (ind2 > FuncLength) RETURN
    
    c=Func(ind2:ind2)
    DO WHILE(c == '*' .OR. c == '/' .OR. c == '%')
      ind2 = CompileUnaryMinus(Comp, Func, ind2+1, Var)

      SELECT CASE(c)
      CASE('*')
        CALL AddCompiledByte(Comp, cMul)
        
      CASE('/')
        CALL AddCompiledByte(Comp, cDiv)
        
      CASE('%')
        CALL AddCompiledByte(Comp, cMod)
        
      END SELECT
      Comp%StackPtr = Comp%StackPtr-1

      IF (ind2 > FuncLength) RETURN
      c=Func(ind2:ind2)
    END DO
  END FUNCTION CompileMult

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileUnaryMinus(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles unary '-'
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    CHARACTER(LEN=1) :: c
    INTEGER :: FuncLength

    FuncLength=LEN_TRIM(Func)

    c = Func(ind:ind)
    IF (c == '-' .OR. c == '!') THEN
      ind2 = ind+1
      IF (ind2 > FuncLength) RETURN
      ind2 = CompilePow(Comp, Func, ind2, Var)

      ! If we are negating a constant, negate the constant itself
      IF (c == '-' .AND. Comp%ByteCode(Comp%ByteCodeSize) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize) = -Comp%Immed(Comp%ImmedSize)
        
        ! If we are negating a negation, we can remove both
      ELSEIF (c == '-' .AND. Comp%ByteCode(Comp%ByteCodeSize) == cNeg) THEN
        CALL RemoveCompiledByte(Comp)
        
      ELSE
        CALL AddCompiledByte(Comp, MERGE(cNeg,cNot,c == '-'))

      END IF
      RETURN
    END IF
    
    ind2 = CompilePow(Comp, Func, ind, Var)
  END FUNCTION CompileUnaryMinus

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompilePow(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles '^'
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    INTEGER :: FuncLength

    FuncLength=LEN_TRIM(Func)

    ind2 = CompileElement(Comp, Func, ind, Var)

    IF (ind2 > FuncLength) RETURN

    DO WHILE(Func(ind2:ind2) == '^')
      ind2 = CompileUnaryMinus(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, cPow)
      Comp%StackPtr = Comp%StackPtr-1
      IF (ind2 > FuncLength) RETURN
    END DO
  END FUNCTION CompilePow

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileElement(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles element
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    CHARACTER(LEN=1) :: c
    REAL(DP) :: rnum
    INTEGER(is) :: n
    INTEGER :: ind1,ind0,ib,in,requiredParams
    LOGICAL :: err

    ind1=ind; c=Func(ind1:ind1)
    IF (c == '(') THEN
      ind1 = CompileExpression(Comp, Func, ind1+1, Var, .FALSE.)
      ind2 = ind1+1   ! Func(ind1:ind1) is ')'
      RETURN
    END IF

    ! Check for numbers
    IF (SCAN(c,'0123456789,') > 0) THEN
      rnum = RealNum (Func(ind1:),ib,in,err)
      IF (err) THEN
        PRINT *, "CompileElement: Invalid number format!"
        CALL sys_halt()
      END IF
      CALL AddImmediate(Comp, rnum)
      CALL AddCompiledByte(Comp, cImmed)
      CALL incStackPtr(Comp)
      ind2 = ind1+in-1
      RETURN
      
    ELSE
      ! Then string must be function, variable or constant
      
      ! Check for mathematical functions
      n = MathFunctionIndex(Func(ind1:))
      IF (n > 0) THEN
        ind2 = ind1+LEN_TRIM(Funcs(n))

        ! Check for IF-THEN-ELSE
        IF (n == cIf) THEN
          ind2 = CompileIf(Comp, Func, ind2+1, Var)
          ! IF-THEN-ELSE cannot be vectorized, note that!
          Comp%isVectorizable = .FALSE.
          RETURN        
        END IF

        requiredParams=MathFunctionParameters(n)        
        ind2 = CompileFunctionParameters(Comp, Func, ind2+1, Var, requiredParams)
        CALL AddFunctionOpcode(Comp, n)
        RETURN
      END IF
      
      ! Check for predefined constant
      n = ConstantIndex(Func(ind1:))
      IF (n > 0) THEN
        ind2 = ind1+LEN_TRIM(Consts(n))+1
        CALL AddImmediate(Comp, Constvals(n))
        CALL AddCompiledByte(Comp, cImmed)
        CALL incStackPtr(Comp)
        RETURN
      END IF

      ! Check for predefined expressions
      n = ExpressionIndex(Func(ind1:))
      IF (n > 0) THEN
        ind2 = ind1+LEN_TRIM(Expressions(n))+1

        ! Recursively compile the given expression
        ind0 = CompileExpression(Comp,ExpressionStrings(n), 1, Var)
        
        ! Proceed with compilation of mathematical function Func afterwards
        RETURN
      END IF

      ! Check for variables
      n = VariableIndex(Func(ind1:), Var, ib, in)
      IF (n > 0) n = VarBegin+n-1
      CALL AddCompiledByte(Comp, n)
      CALL incStackPtr(Comp)
      ind2 = ind1+in-1
      RETURN
    END IF

    PRINT *, "CompileElement: An unexpected error occured."
    CALL sys_halt()
  END FUNCTION CompileElement

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileFunctionParameters(Comp, Func, ind, Var,&
      & requiredParams) RESULT(ind2)

!<description>
    ! Compiles function parameters
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var

    ! Number of required parameters
    INTEGER, INTENT(IN) :: requiredParams
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    INTEGER :: curStackPtr

    ind2 = ind
    IF (requiredParams > 0) THEN
      
      curStackPtr = Comp%StackPtr
      ind2 = CompileExpression(Comp, Func, ind, Var, .FALSE.)

      IF (Comp%StackPtr /= curStackPtr+requiredParams) THEN
        PRINT *, "CompileFunctionParameters: Illegal number of parameters to function!"
        CALL sys_halt()
      END IF

      Comp%StackPtr = Comp%StackPtr-(requiredParams-1)

    ELSE
      
      CALL incStackPtr(Comp)
      
    END IF
    ind2=ind2+1
  END FUNCTION CompileFunctionParameters

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileIf(Comp, Func, ind, Var) RESULT(ind2)

!<description>
    ! Compiles if()
!</description>

!<input>
    ! Function substring
    CHARACTER(LEN=*), INTENT(IN) :: Func

    ! Begin position substring
    INTEGER, INTENT(IN) :: ind

    ! Array with variable names
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE(t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    INTEGER :: ind2
!</result>
!</function>

    ! local variables
    INTEGER :: FuncLength,curByteCodeSize,curByteCodeSize2,curImmedSize2

    FuncLength=LEN_TRIM(Func)

    ind2 = CompileExpression(Comp, Func, ind, Var, .TRUE.) ! Condition branch
    IF (ind2 > FuncLength) RETURN

    IF (Func(ind2:ind2) /= ',') THEN
      PRINT *, "CompileIf: Illegal number of parameters to function!"
      CALL sys_halt()
    END IF
    CALL AddCompiledByte(Comp, cIf)
    curByteCodeSize = Comp%ByteCodeSize
    CALL AddCompiledByte(Comp, 0_is) ! Jump index will be set below
    CALL AddCompiledByte(Comp, 0_is) ! Immed jump index will be set below
    Comp%StackPtr = Comp%StackPtr-1

    ind2 = CompileExpression(Comp, Func, ind2+1, Var, .TRUE.) ! Then branch
    IF (ind2 > FuncLength) RETURN

    IF (Func(ind2:ind2) /= ',') THEN
      PRINT *, "CompileIf: Illegal number of parameters to function!"
      CALL sys_halt()
    END IF
    CALL AddCompiledByte(Comp, cJump)
    curByteCodeSize2 = Comp%ByteCodeSize
    curImmedSize2 = Comp%ImmedSize
    CALL AddCompiledByte(Comp, 0_is) ! Jump index will be set below
    CALL AddCompiledByte(Comp, 0_is) ! Immed jump index will be set below
    Comp%StackPtr = Comp%StackPtr-1

    ind2 = CompileExpression(Comp, Func, ind2+1, Var, .TRUE.) ! Else branch
    IF (ind2 > FuncLength) RETURN

    IF (Func(ind2:ind2) /= ')') THEN
      PRINT *, "CompileIf: Illegal number of parameters to function!"
      CALL sys_halt()
    END IF
    
    ! Set jump indices
    IF (ASSOCIATED(Comp%ByteCode)) THEN
      Comp%ByteCode(curByteCodeSize+1)  = curByteCodeSize2+2
      Comp%ByteCode(curByteCodeSize+2)  = curImmedSize2+1
      Comp%ByteCode(curByteCodeSize2+1) = Comp%ByteCodeSize
      Comp%ByteCode(curByteCodeSize2+2) = Comp%ImmedSize+1
    END IF
    
    ind2=ind2+1
  END FUNCTION CompileIf

  ! *****************************************************************************

!<function>

  ELEMENTAL FUNCTION DbleTOLogc(d) RESULT(l)

!<description>
    ! This function transforms a Double into a Logical
!</description>

!<input>
    ! Double variable
    REAL(DP), INTENT(IN) :: d
!</input>

!<result>
    ! Logical variable
    LOGICAL :: l
!</result>
!</function>
    
    l=MERGE(.TRUE.,.FALSE.,ABS(1-d)<=1e-12)
  END FUNCTION DbleTOLogc

  ! *****************************************************************************

!<function>
  
  ELEMENTAL FUNCTION LogcToDble(l) RESULT(d)

!<description>
    ! This function transforms a Logical into a Double
!</description>

!<input>
    ! Logical variable
    LOGICAL, INTENT(IN) :: l
!</input>

!<result>
    ! Double variable
    REAL(DP) :: d
!</result>
!</function>
      
    d=MERGE(1._DP,0._DP,l)
  END FUNCTION LogcToDble

  ! *****************************************************************************

!<function>
  
  ELEMENTAL FUNCTION DegToRad(d) RESULT(r)

!<description>
    ! This function converts DEG to RAD
!</description>

!<input>
    ! DEG
    REAL(DP), INTENT(IN) :: d
!</input>

!<result>
    ! RAD
    REAL(DP) :: r
!</result>
!</function>
      
    r=d*(3.141592653589793115997963468544185161590576171875_DP/180._DP)
  END FUNCTION DegToRad

  ! *****************************************************************************

!<function>

  ELEMENTAL FUNCTION RadToDeg(r) RESULT(d)
    
!<description>
    ! This function converts RAD to DEG
!</description>

!<input>
    ! RAD
    REAL(DP), INTENT(IN) :: r
!</input>

!<result>
    ! DEG
    REAL(DP) :: d
!</result>
!</function>
    
    d=r*(180._DP/3.141592653589793115997963468544185161590576171875_DP)
  END FUNCTION RadToDeg
END MODULE fparser
