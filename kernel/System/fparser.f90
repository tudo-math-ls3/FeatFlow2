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
!# The following routines can be found in this module:
!#
!# 1.) fparser_create
!#     -> Create function parser
!#
!# 2.) fparser_release
!#     -> Release function parser
!#
!# 3.) fparser_parseFunction
!#     -> Parse function string and compile ot into bytecode
!#
!# 4.) fparser_evalFunction = fparser_evalFunctionScalar /
!#                            fparser_evalFunctionArray
!#     -> Evaluate bytecode
!#
!# 5.) fparser_ErrorMsg
!#     -> Get error message from function parser
!#
!# 6.) fparser_PrintByteCode
!#     -> Print the bytecode stack (very technical!)
!#
!# The following internal routines can be found in this module:
!#
!# 1.) CheckSyntax
!#     -> Check syntax of function string before compiling bytecode
!#
!# 2.) ParseErrMsg
!#     -> Print detailed error message and terminate
!#
!# 3.) isOperator
!#     -> Return size of operator and 0 otherwise
!#
!# 4.) MathFunctionIndex
!#     -> Return index of mathematical function and 0 otherwise
!#
!# 5.) MathFunctionParameters
!#     -> Return number of required function parameters
!#
!# 6.) ConstantIndex
!#     -> Return index of predefined constant and 0 otherwise
!#
!# 7.) VariableIndex
!#     -> Return index of variable
!#
!# 8.) RemoveSpaces
!#     -> Remove spaces from string
!#
!# 8.) Replace
!#     -> Replace all appearances of one character set by another
!#        character set in a given string
!#
!# 9.) Compile
!#     -> Compile function string into bytecode
!#
!# 10.) incStackPtr
!#     -> Increase stack pointer
!#
!# 11.) AddCompiledByte
!#     -> Add compiled byte to bytecode stack
!#
!# 11.) RemoveCompiledByte
!#     -> Remove last compiled byte from bytecode stack
!#
!# 12.) AddImmediate
!#     -> Add immediate to immediate stack
!#
!# 13.) AddFunctionOpcode
!#     -> Add function opcode to bytecode stack
!#
!# 14.) RealNum
!#     -> Get real number from string
!#
!# 15.) LowCase
!#     -> Transform upper case letters to lower case letters
!#
!# 16.) CompileExpression
!#      -> Compile ','
!#
!# 17.) CompileOr
!#      -> Compile '|'
!#
!# 18.) CompileAnd
!#      -> Compile '&'
!#
!# 19.) CompileComparison
!#      -> Compile '=', '<', and '>'
!#
!# 20.) CompileAddition
!#      -> Compile '+' and '-'
!#
!# 21.) CompileMult
!#      -> Compile '*', '/', and '%'
!#
!# 22.) CompileUnaryMinus
!#      -> Compile unary '-'
!#
!# 23.) CompilePow
!#      -> Compile '^'
!#
!# 24.) CompileElement
!#      -> Compile mathematical function, variable, constant and number
!#
!# 25.) CompileFunctionParameters
!#      -> Compile function parameters
!#
!# 26.) CompileIf
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
  PUBLIC :: fparser_create
  PUBLIC :: fparser_release
  PUBLIC :: fparser_parseFunction
  PUBLIC :: fparser_evalFunction
  PUBLIC :: fparser_ErrorMsg
  PUBLIC :: fparser_PrintByteCode
  public :: lowcase

  INTERFACE fparser_evalFunction
    MODULE PROCEDURE fparser_evalFunctionScalar
    MODULE PROCEDURE fparser_evalFunctionArray
  END INTERFACE

!<constants>

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
                            cAsin        = 39, &
                            cAcos        = 40
  INTEGER(is), PARAMETER :: cAtan        = 41, &
                            cAcosh       = 42, &
                            cAsinh       = 43, &
                            cAtanh       = 44, &
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
                                                                  'asin ', &
                                                                  'acos ', &
                                                                  'atan ', &
                                                                  'acosh', &
                                                                  'acinh', &
                                                                  'atanh', &
                                                                  'ceil ', &
                                                                  'floor', &
                                                                  'csc  ', &
                                                                  'sec  ' /)
!</constantblock>

!<constantblock description="constants for parser">
  INTEGER(is), PARAMETER :: C_PI   = 1, &
                            C_EXP  = 2
!</constantblock>

!<constantblock description="constant names for parser">
  CHARACTER(LEN=3), DIMENSION(C_PI:C_EXP), PARAMETER :: Consts = (/ 'pi ', &
                                                                    'exp' /)
!</constantblock>

!<constantblock description="constant values for parser">
  REAL(DP), DIMENSION(C_PI:c_EXP), PARAMETER :: Constvals = (/&
      & 3.141592653589793115997963468544185161590576171875_DP, &
      & 2.718281828459045090795598298427648842334747314453125_DP /)
!</constantblock>

!</constants>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<globals>

  ! Associates function strings
  INTEGER, PRIVATE, DIMENSION(:),  POINTER :: ipos => NULL()
!</globals>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<types>

!<typeblock>
  
  ! Type block for storing all information of the function parser
  TYPE t_fparser
    PRIVATE
    
    ! Array of function parser components
    TYPE(t_fparserComponent), DIMENSION(:), POINTER :: Comp => NULL()
    
    ! Number of parser components
    INTEGER :: nComp = 0
    
    ! Evaluation error code: =0: no error occured, >0
    INTEGER :: EvalErrType = 0
  END TYPE t_fparser

!</typeblock>

!<typeblock>

  ! Type block for storing the bytecode of the function parser for one component
  TYPE t_fparserComponent
    ! Size of bytecode
    INTEGER :: ByteCodeSize = 0

    ! Size of immediates
    INTEGER :: ImmedSize = 0

    ! Stack size
    INTEGER :: StackSize = 0

    ! Stack pointer
    INTEGER :: StackPtr = 0

    ! Use degree conversion for some functions
    LOGICAL :: useDegreeConversion = .FALSE.
    
    ! Bytecode
    INTEGER(is), DIMENSION(:), POINTER :: ByteCode => NULL()
    
    ! Immediates
    REAL(DP), DIMENSION(:), POINTER :: Immed => NULL()
  END TYPE t_fparserComponent
!</typeblock>

!</types>


CONTAINS
  
  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE fparser_create (rparser,n)

!<description>
    ! Initialize function parser for n functions
!</description>

!<input>
    ! Number of functions
    INTEGER, INTENT(IN) :: n                                 
!</input>

!<output>
    ! Function parser object
    TYPE(t_fparser) :: rparser
!</output>
!</subroutine>
    
    ALLOCATE (rparser%Comp(n))
    rparser%nComp=n
  END SUBROUTINE fparser_create

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_release (rparser)

!<description>
    ! Release function parser
!</description>

!<inputoutput>
    ! Function parser
    TYPE(t_fparser), INTENT(INOUT) :: rparser
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER :: i
    
    DO i=1,rparser%nComp
      IF (ASSOCIATED(rparser%Comp(i)%ByteCode)) DEALLOCATE(rparser%Comp(i)%ByteCode)
      IF (ASSOCIATED(rparser%Comp(i)%Immed)) DEALLOCATE(rparser%Comp(i)%Immed)
    END DO
    DEALLOCATE(rparser%Comp)
    rparser%nComp=0
  END SUBROUTINE fparser_release
  
  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE fparser_parseFunction (rparser, i, FuncStr, Var, useDegrees)

!<description>
    ! Parse ith function string FuncStr and compile it into bytecode
!</description>
    
!<input>
    ! Function identifier
    INTEGER, INTENT(IN) :: i

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
    
    IF (i < 1 .OR. i > rparser%nComp) THEN
      WRITE(*,*) '*** Parser error: Function number ',i,' out of range'
      STOP
    END IF
    
    ! Char. positions in orig. string
    ALLOCATE (ipos(LEN_TRIM(FuncStr)))

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

    ! Check for valid syntax
    CALL CheckSyntax (Func,FuncStr,Var)
    
    ! Compile into bytecode
    IF (PRESENT(useDegrees)) rparser%Comp(i)%useDegreeConversion&
        &=useDegrees
    CALL Compile (rparser%Comp(i),Func,Var)

    DEALLOCATE (ipos)
  END SUBROUTINE fparser_parseFunction

  ! *****************************************************************************

!<function>

  FUNCTION fparser_evalFunctionScalar (rparser, i, Val, EvalErrType) RESULT (res)

!<description>
    ! Evaluate bytecode of ith function for the values passed in
    ! array Val(:)
!</description>

!<input>
    ! Function identifier
    INTEGER, INTENT(IN) :: i

    ! Variable values
    REAL(DP), DIMENSION(:), INTENT(IN) :: Val
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_fparser),  INTENT(INOUT) :: rparser
!</inputoutput>

!<output>
    ! OPTIONAL Error code for function evaluation
    INTEGER, INTENT(OUT), OPTIONAL :: EvalErrType
!</output>

!<result>
    ! Evaluated function
    REAL(DP) :: res
!</result>
!</function>

    ! local parameters
    REAL(DP), PARAMETER :: zero = 0._DP
    
    ! local variables
    TYPE(t_fparserComponent), POINTER :: Comp
    INTEGER :: IP, DataPtr, StackPtr,jumpAddr,immedAddr,h_Stack
    REAL(DP), DIMENSION(:), POINTER :: Stack
    REAL(DP) :: daux
    
    ! Initialization
    DataPtr  = 1
    StackPtr = 0
    IP       = 0
    Comp => rparser%Comp(i)

    ! Allocate stack
    CALL storage_new('fparser_evalFunctionScalar','Stack',INT(Comp%StackSize+1,I32),&
        ST_DOUBLE,h_Stack,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_double(h_Stack,Stack)

    DO WHILE(IP < Comp%ByteCodeSize)
      IP=IP+1
      
      SELECT CASE (Comp%ByteCode(IP))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      CASE (cAbs)
        Stack(StackPtr)=ABS(Stack(StackPtr))
        
      CASE (cAcos)
        IF ((Stack(StackPtr) < -1._DP) .OR. &
            (Stack(StackPtr) > 1._DP)) THEN
          rparser%EvalErrType=4; res=zero
          CALL storage_free(h_Stack)
          RETURN
        ENDIF
        Stack(StackPtr)=ACOS(Stack(StackPtr))

      CASE (cAsin)
        IF ((Stack(StackPtr) < -1._DP) .OR. &
            (Stack(StackPtr) > 1._DP)) THEN
          rparser%EvalErrType=4; res=zero
          CALL storage_free(h_Stack)
          RETURN
        ENDIF
        Stack(StackPtr)=ASIN(Stack(StackPtr))

      CASE (cAtan)
        Stack(StackPtr)=ATAN(Stack(StackPtr))

      CASE  (cAtan2)
        Stack(StackPtr-1)=ATAN2(Stack(StackPtr -1),Stack(StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cAcosh)
        daux=Stack(StackPtr)+SQRT(Stack(StackPtr)**2-1)
        IF (daux <= 0._DP) THEN
          rparser%EvalErrType=5; res=zero
          CALL storage_free(h_Stack)
          RETURN
        END IF
        Stack(StackPtr)=LOG(daux)
        
      CASE (cAnint)
        Stack(StackPtr)=ANINT(Stack(StackPtr))

      CASE (cAint)
        Stack(StackPtr)=AINT(Stack(StackPtr))

      CASE (cAsinh)
        daux=Stack(StackPtr)+SQRT(Stack(StackPtr)+1)
        IF (daux <= 0._DP) THEN
          rparser%EvalErrType=5; res=zero
          CALL storage_free(h_Stack)
          RETURN
        END IF
        Stack(StackPtr)=LOG(daux)

      CASE (cAtanh)
        IF (Stack(StackPtr) == -1._DP) THEN
          rparser%EvalErrType=6; res=zero
          CALL storage_free(h_Stack)
          RETURN
        END IF
        daux=(1+Stack(StackPtr))/(1-Stack(StackPtr))
        IF (daux <= 0._DP) THEN
          rparser%EvalErrType=3; res=zero
          CALL storage_free(h_Stack)
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
        IF (daux == 0) THEN 
          rparser%EvalErrType=1; res=zero
          CALL storage_free(h_Stack)
          RETURN
        END IF
        Stack(StackPtr)=1._DP/daux

      CASE (cCsc)
        daux=SIN(Stack(StackPtr))
        IF (daux==0._DP) THEN
          rparser%EvalErrType=1; res=zero
          CALL storage_free(h_Stack)
          RETURN
        ENDIF
        Stack(StackPtr)=1._DP/daux

      CASE (cExp)
        Stack(StackPtr)=EXP(Stack(StackPtr))

      CASE (cFloor)
        Stack(StackPtr)=FLOOR(Stack(StackPtr))

      CASE (cIf)
        IP=IP+1; jumpAddr = Comp%ByteCode(IP)
        IP=IP+1; immedAddr = Comp%ByteCode(IP)
        IF (.NOT.DbleToLogc(Stack(StackPtr))) THEN
          IP      = jumpAddr
          DataPtr = immedAddr
        END IF
        StackPtr=StackPtr-1
        
      CASE (cLog)
        IF (Stack(StackPtr) <= 0._DP) THEN
          rparser%EvalErrType=3; res=zero
          CALL storage_free(h_Stack)
          RETURN
        ENDIF
        Stack(StackPtr)=LOG(Stack(StackPtr)) 

      CASE (cLog10)
        IF (Stack(StackPtr) <= 0._DP) THEN
          rparser%EvalErrType=3; res=zero
          CALL storage_free(h_Stack)
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
        IF (daux==0._DP) THEN
          rparser%EvalErrType=1; res=zero
          CALL storage_free(h_Stack)
          RETURN
        ENDIF
        Stack(StackPtr)=1._DP/daux
        
      CASE (cSin)
        Stack(StackPtr)=SIN(Stack(StackPtr))
        
      CASE(cSinh)
        Stack(StackPtr)=SINH(Stack(StackPtr))
      
      CASE(cSqrt)
        IF (Stack(StackPtr) < 0._DP) THEN
          rparser%EvalErrType=3; res=zero
          CALL storage_free(h_Stack)
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
        StackPtr=StackPtr+1
        Stack(StackPtr)=Comp%Immed(DataPtr)
        DataPtr=DataPtr+1

      CASE (cJump)
        DataPtr=Comp%ByteCode(IP+2)
        IP     =Comp%ByteCode(IP+1)

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
        IF (Stack(StackPtr)==0._DP) THEN
          rparser%EvalErrType=1; res=zero
          CALL storage_free(h_Stack)
          RETURN
        ENDIF
        Stack(StackPtr-1)=Stack(StackPtr-1)/Stack(StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cMod)
        IF (Stack(StackPtr)==0._DP) THEN
          rparser%EvalErrType=1; res=zero
          CALL storage_free(h_Stack)
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
        Stack(StackPtr)=Val(Comp%ByteCode(IP)-VarBegin+1)
      END SELECT
    END DO
    rparser%EvalErrType = 0
    res = Stack(StackPtr)

    ! Set error code (if required)
    IF (PRESENT(EvalErrType)) EvalErrType = rparser%EvalErrType

    ! Free memory
    CALL storage_free(h_Stack)
  END FUNCTION fparser_evalFunctionScalar

  ! *****************************************************************************

!<function>

  FUNCTION fparser_evalFunctionArray (rparser, i, isize, Ishape, Val, EvalErrType) RESULT (Res)

!<description>
    ! Evaluate bytecode of ith function for an array of values passed
    ! in Val(:,:)
!</description>

!<input>
    ! Function identifier
    INTEGER, INTENT(IN) :: i

    ! Size of the resulting array
    INTEGER, INTENT(IN) :: isize

    ! Shape of the value array
    INTEGER, DIMENSION(2), INTENT(IN) :: Ishape

    ! Variable values
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Val
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_fparser),  INTENT(INOUT) :: rparser
!</inputoutput>

!<output>
    ! OPTIONAL Error code for function evaluation
    INTEGER, INTENT(OUT), OPTIONAL :: EvalErrType
!</output>

!<result>
    ! Evaluated function
    REAL(DP), DIMENSION(isize) :: Res
!</result>
!</function>

    ! local parameters
    REAL(DP), PARAMETER :: zero = 0._DP
    
    ! local variables
    TYPE(t_fparserComponent), POINTER :: Comp
    INTEGER :: IP,DataPtr,StackPtr,h_Stack,j
    REAL(DP), DIMENSION(:,:), POINTER :: Stack
    REAL(DP) :: daux
    LOGICAL :: bfirstdim
    INTEGER(I32), DIMENSION(2) :: IsizeAlloc

    ! Check if dimensions agree
    IF (ANY(SHAPE(Val) /= Ishape)) THEN
      PRINT *, 'fparser_evalFunctionArray: invalid dimensions!'
      STOP
    END IF
    
    IF (isize == Ishape(1)) THEN
      bfirstdim=.TRUE.
    ELSEIF (isize == Ishape(2)) THEN
      bfirstdim=.FALSE.
    ELSE
      PRINT *, 'fparser_evalFunctionArray: invalid dimensions!'
      STOP
    END IF
    
    ! Initialization
    DataPtr  = 1
    StackPtr = 0
    IP       = 0
    Comp => rparser%Comp(i)
    rparser%EvalErrType=0

    ! Allocate stack
    IsizeAlloc = (/INT(isize,I32),INT(Comp%StackSize+1,I32)/)
    CALL storage_new('fparser_evalFunctionArray','Stack',&
        IsizeAlloc,ST_DOUBLE,h_Stack,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_double2D(h_Stack,Stack)

    DO WHILE(IP < Comp%ByteCodeSize)
      IP=IP+1
      
      SELECT CASE (Comp%ByteCode(IP))
        !------------------------------------------------------------
        ! Functions
        !------------------------------------------------------------
      CASE (cAbs)
        Stack(:,StackPtr)=ABS(Stack(:,StackPtr))
        
      CASE (cAcos)
!$omp parallel do default(shared) private(j)
        DO j=1,isize
          IF ((Stack(j,StackPtr) < -1._DP) .OR.&
              (Stack(j,StackPtr) > 1._DP)) THEN
            rparser%EvalErrType=4; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=ACOS(Stack(j,StackPtr))
          END IF
        END DO
!$omp end parallel do

      CASE (cAsin)
!$omp parallel do default(shared) private(j)
        DO j=1,isize
          IF ((Stack(j,StackPtr) < -1._DP) .OR.&
              (Stack(j,StackPtr) > 1._DP)) THEN
            rparser%EvalErrType=4; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=ASIN(Stack(j,StackPtr))
          END IF
        END DO
!$omp end parallel do
        
      CASE (cAtan)
        Stack(:,StackPtr)=ATAN(Stack(:,StackPtr))

      CASE  (cAtan2)
        Stack(:,StackPtr-1)=ATAN2(Stack(:,StackPtr -1),Stack(:,StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cAcosh)
!$omp parallel do default(shared) private(j,daux)
        DO j=1,isize
          daux=Stack(j,StackPtr)+SQRT(Stack(j,StackPtr)**2-1)
          IF (daux <= 0._DP) THEN
            rparser%EvalErrType=5; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=LOG(daux)
          END IF
        END DO
!$omp end parallel do
        
      CASE (cAnint)
        Stack(:,StackPtr)=ANINT(Stack(:,StackPtr))

      CASE (cAint)
        Stack(:,StackPtr)=AINT(Stack(:,StackPtr))

      CASE (cAsinh)
!$omp parallel do default(shared) private(j,daux)
        DO j=1,isize
          daux=Stack(j,StackPtr)+SQRT(Stack(j,StackPtr)+1)
          IF (daux <= 0._DP) THEN
            rparser%EvalErrType=5; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=LOG(daux)
          END IF
        END DO
!$omp end parallel do

      CASE (cAtanh)
!$omp parallel do default(shared) private(j,daux)
        DO j=1,isize
          IF (Stack(j,StackPtr) == -1._DP) THEN
            rparser%EvalErrType=6; Res(j)=zero
          END IF
          daux=(1+Stack(j,StackPtr))/(1-Stack(j,StackPtr))
          IF (daux <= 0._DP) THEN
            rparser%EvalErrType=3; Res(j)=zero
          ELSE
            Stack(j,StackPtr) = LOG(daux)/2._DP
          END IF
        END DO
!$omp end parallel do
        
      CASE (cCeil)
        Stack(:,StackPtr)=CEILING(Stack(:,StackPtr))

      CASE (cCos)
        Stack(:,StackPtr)=COS(Stack(:,StackPtr)) 
        
      CASE (cCosh)
        Stack(:,StackPtr)=COSH(Stack(:,StackPtr))

      CASE (cCot)
!$omp parallel do default(shared) private(j,daux)
        DO j=1,isize
          daux=TAN(Stack(j,StackPtr))
          IF (daux == 0) THEN 
            rparser%EvalErrType=1; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=1._DP/daux
          END IF
        END DO
!$omp end parallel do

      CASE (cCsc)
!$omp parallel do default(shared) private(j,daux)
        DO j=1,isize
          daux=SIN(Stack(j,StackPtr))
          IF (daux==0._DP) THEN
            rparser%EvalErrType=1; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=1._DP/daux
          END IF
        END DO
!$omp end parallel do

      CASE (cExp)
        Stack(:,StackPtr)=EXP(Stack(:,StackPtr))

      CASE (cFloor)
        Stack(:,StackPtr)=FLOOR(Stack(:,StackPtr))

      CASE (cIf)
        ! At the moment, IF-THEN-ELSE cannot be handled for multiple
        ! values. Hence, stop this subroutine an process each item separately
        CALL storage_free(h_Stack)
        IF (bfirstdim) THEN
          DO j=1,isize
            Res(j)=fparser_evalFunction(rparser,i,Val(j,:),EvalErrType)
          END DO
        ELSE
          DO j=1,isize
            Res(j)=fparser_evalFunction(rparser,i,Val(:,j),EvalErrType)
          END DO
        END IF
        CALL storage_free(h_Stack)
        RETURN
        
      CASE (cLog)
!$omp parallel do default(shared) private(j)
        DO j=1,isize
          IF (Stack(j,StackPtr) <= 0._DP) THEN
            rparser%EvalErrType=3; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=LOG(Stack(j,StackPtr)) 
          END IF
        END DO
!$omp end parallel do

      CASE (cLog10)
!$omp parallel do default(shared) private(j)
        DO j=1,isize
          IF (Stack(j,StackPtr) <= 0._DP) THEN
            rparser%EvalErrType=3; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=LOG10(Stack(j,StackPtr))
          END IF
        END DO
!$omp end parallel do

      CASE (cMax)
        Stack(:,StackPtr-1)=MAX(Stack(:,StackPtr-1),Stack(:,StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cMin)
        Stack(:,StackPtr-1)=MIN(Stack(:,StackPtr-1),Stack(:,StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cSec)
!$omp parallel do default(shared) private(j,daux)
        DO j=1,isize
          daux=COS(Stack(j,StackPtr))
          IF (daux==0._DP) THEN
            rparser%EvalErrType=1; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=1._DP/daux
          END IF
        END DO
!$omp end parallel do
        
      CASE (cSin)
        Stack(:,StackPtr)=SIN(Stack(:,StackPtr))
        
      CASE(cSinh)
        Stack(:,StackPtr)=SINH(Stack(:,StackPtr))
      
      CASE(cSqrt)
!$omp parallel do default(shared) private(j)
        DO j=1,isize
          IF (Stack(j,StackPtr) < 0._DP) THEN
            rparser%EvalErrType=3; Res(j)=zero
          ELSE
            Stack(j,StackPtr)=SQRT(Stack(j,StackPtr))
          END IF
        END DO
!$omp end parallel do

      CASE (cTan)
        Stack(:,StackPtr)=TAN(Stack(:,StackPtr))
        
      CASE (cTanh)
        Stack(:,StackPtr)=TANH(Stack(:,StackPtr))
        
        !------------------------------------------------------------
        ! Misc
        !------------------------------------------------------------
      CASE (cImmed)
        StackPtr=StackPtr+1
        Stack(:,StackPtr)=Comp%Immed(DataPtr)
        DataPtr=DataPtr+1

      CASE (cJump)
        DataPtr=Comp%ByteCode(IP+2)
        IP     =Comp%ByteCode(IP+1)

        !------------------------------------------------------------
        ! Operators
        !------------------------------------------------------------
      CASE (cNeg)
        Stack(:,StackPtr)=-Stack(:,StackPtr)
        
      CASE (cAdd)
        Stack(:,StackPtr-1)=Stack(:,StackPtr-1)+Stack(:,StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cSub)
        Stack(:,StackPtr-1)=Stack(:,StackPtr-1)-Stack(:,StackPtr)
        StackPtr=StackPtr-1

      CASE (cMul)
        Stack(:,StackPtr-1)=Stack(:,StackPtr-1)*Stack(:,StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cDiv)
!$omp parallel do default(shared) private(j)
        DO j=1,isize
          IF (Stack(j,StackPtr)==0._DP) THEN
            rparser%EvalErrType=1; Res(j)=zero
          ELSE
            Stack(j,StackPtr-1)=Stack(j,StackPtr-1)/Stack(j,StackPtr)
          END IF
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1
        
      CASE (cMod)
!$omp parallel do default(shared) private(j)
        DO j=1,isize
          IF (Stack(j,StackPtr)==0._DP) THEN
            rparser%EvalErrType=1; Res(j)=zero
          ELSE
            Stack(j,StackPtr-1)=MOD(Stack(j,StackPtr-1),Stack(j,StackPtr))
          END IF
        END DO
!$omp end parallel do
        StackPtr=StackPtr-1

      CASE (cPow)
        Stack(:,StackPtr-1)=Stack(:,StackPtr-1)**Stack(:,StackPtr)
        StackPtr=StackPtr-1
        
      CASE (cEqual)
        Stack(:,StackPtr-1)=LogcToDble(Stack(:,StackPtr-1) == Stack(:,StackPtr))
        StackPtr=StackPtr-1

      CASE (cNEqual)
        Stack(:,StackPtr-1)=LogcToDble(Stack(:,StackPtr-1) /= Stack(:,StackPtr))
        StackPtr=StackPtr-1

      CASE (cLess)
        Stack(:,StackPtr-1)=LogcToDble(Stack(:,StackPtr-1) < Stack(:,StackPtr))
        StackPtr=StackPtr-1

      CASE (cLessOrEq)
        Stack(:,StackPtr-1)=LogcToDble(Stack(:,StackPtr-1) <= Stack(:,StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cGreater)
        Stack(:,StackPtr-1)=LogcToDble(Stack(:,StackPtr-1) > Stack(:,StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cGreaterOrEq)
        Stack(:,StackPtr-1)=LogcToDble(Stack(:,StackPtr-1) >= Stack(:,StackPtr))
        StackPtr=StackPtr-1
        
      CASE (cAnd)
        Stack(:,StackPtr-1)=LogcToDble(DbleToLogc(Stack(:,StackPtr-1)) .AND. &
            DbleToLogc(Stack(:,StackPtr)) )
        StackPtr=StackPtr-1

      CASE (cOr)
        Stack(:,StackPtr-1)=LogcToDble(DbleToLogc(Stack(:,StackPtr-1)) .OR. &
            DbleToLogc(Stack(:,StackPtr)) )
        StackPtr=StackPtr-1

      CASE (cNot)
        Stack(:,StackPtr)=LogcToDble( .NOT. DbleToLogc(Stack(:,StackPtr)) )
        
        !------------------------------------------------------------
        ! Degrees-radians conversion
        !------------------------------------------------------------
      CASE (cDeg)
        Stack(:,StackPtr)=RadToDeg(Stack(:,StackPtr))
        
      CASE (cRad)
        Stack(:,StackPtr)=DegToRad(Stack(:,StackPtr))
        
      CASE DEFAULT
        StackPtr=StackPtr+1
        IF (bfirstdim) THEN
          Stack(:,StackPtr)=Val(:,Comp%ByteCode(IP)-VarBegin+1)
        ELSE
          Stack(:,StackPtr)=Val(Comp%ByteCode(IP)-VarBegin+1,:)
        END IF
      END SELECT

      ! Check if an error occured
      IF (rparser%EvalErrType /= 0) THEN
        CALL storage_free(h_Stack)    
        RETURN
      END IF

    END DO
    rparser%EvalErrType = 0
    Res = Stack(:,StackPtr)

    ! Set error code (if required)
    IF (PRESENT(EvalErrType)) EvalErrType = rparser%EvalErrType

    ! Free memory
    CALL storage_free(h_Stack)    
  END FUNCTION fparser_evalFunctionArray  
  
  ! *****************************************************************************

!<function>

  FUNCTION fparser_ErrorMsg (rparser) RESULT (msg)

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
    ! Function parser
    TYPE(t_fparser), INTENT(IN) :: rparser
!</input>   
    
!<result>
    ! Error messages
    CHARACTER (LEN=LEN(m)) :: msg
!</result>
!</function>

    
    
    IF (rparser%EvalErrType < 1 .OR. rparser%EvalErrType > SIZE(m)) THEN
      msg = ''
    ELSE
      msg = m(rparser%EvalErrType)
    ENDIF
  END FUNCTION fparser_ErrorMsg

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_PrintByteCode(rparser, i)

!<description>
    ! Print the compiled bytecode stack
!</description>

!<input>
    ! Function parser
    TYPE(t_fparser), INTENT(IN) :: rparser

    ! Function identifier
    INTEGER, INTENT(IN) :: i
!</input>
!</subroutine>

    ! local variables
    TYPE(t_fparserComponent), POINTER :: Comp
    CHARACTER(LEN=5) :: n
    INTEGER :: IP, DataPtr, StackPtr,nparams
    INTEGER(is) :: opCode
    
    nparams=1
    DataPtr  = 1
    StackPtr = 0
    IP       = 0
    Comp => rparser%Comp(i)

    DO WHILE(IP < Comp%ByteCodeSize)
      IP=IP+1
      
      WRITE(*,FMT='(I8.8,1X,":",1X)',ADVANCE="NO") IP
      
      opCode=Comp%ByteCode(IP)
      SELECT CASE(opCode)
        
      CASE(cIf)
        WRITE(*,FMT='(A,1X,T10,I8.8)') "jz",Comp%ByteCode(IP+1)+1
        IP=IP+2

      CASE(cJump)
        WRITE(*,FMT='(A,1X,T10,I8.8)') "jump",Comp%ByteCode(IP+1)+1
        IP=IP+2

      CASE(cImmed)
        WRITE(*,FMT='(A,1X,T10,ES8.2)') "push",Comp%Immed(DataPtr)
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
            n = Funcs(opCode)
            nparams=MathFunctionParameters(opCode)
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
  
  SUBROUTINE CheckSyntax (Func,FuncStr,Var)

!<description>
    ! Check syntax of function string, returns 0 if syntax is ok
!</description>
    
!<input>
    ! Function string without spaces
    CHARACTER (LEN=*), INTENT(in) :: Func

    ! Original function string
    CHARACTER (LEN=*), INTENT(in) :: FuncStr

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
    INTEGER :: ParCnt,Ind,Ind2,ib,in,lFunc,idummy
    
    Ind = 1
    ParCnt = 0
    lFunc = LEN_TRIM(Func)
    CALL stack_create(functionParenthDepth,MAX(5,INT(lFunc/4._DP)),ST_INT)

    DO
      IF (Ind > lFunc) CALL ParseErrMsg (Ind, FuncStr)
      c = Func(Ind:Ind)

      ! Check for valid operand (must appear)

      ! Check for leading - or !
      IF (c == '-' .OR. c == '!') THEN                      
        Ind = Ind+1
        IF (Ind > lFunc) CALL ParseErrMsg (Ind-1, FuncStr, 'Premature &
            &end of string')
        c = Func(Ind:Ind)
      END IF
            
      ! Check for math function
      n = MathFunctionIndex (Func(Ind:))
      IF (n > 0) THEN
        ! Math function found
        Ind = Ind+LEN_TRIM(Funcs(n))
        IF (Ind > lFunc) CALL ParseErrMsg (Ind-1, FuncStr, 'Premature &
            &end of string')
        c = Func(Ind:Ind)
        IF (c /= '(') CALL ParseErrMsg (Ind, FuncStr, 'Expecting ( aft&
            &er function')

        Ind2=Ind+1
        IF (Ind2 > lFunc) CALL ParseErrMsg (Ind, FuncStr, 'Premature e&
            &nd of string')
        IF (Func(Ind2:Ind2) == ')') THEN
          Ind=Ind2+1
          IF (Ind > lFunc) CALL ParseErrMsg (Ind2, FuncStr, 'Premature&
              & end of string')
          c = Func(Ind:Ind)
          
          ! Ugly, but other methods would just be uglier ...
          GOTO 999
        END IF
        
        ! Push counter for parenthesss to stack
        CALL stack_pushback(functionParenthDepth,ParCnt+1)
      END IF
      
      ! Check for opening parenthesis
      IF (c == '(') THEN
        ParCnt = ParCnt+1
        Ind = Ind+1
        IF (Ind > lFunc) CALL ParseErrMsg (Ind-1, FuncStr, 'Premature &
            &end of string')
        IF (Func(Ind:Ind) == ')') CALL ParseErrMsg (Ind, FuncStr, 'Emp&
            &ty parantheses')
        CYCLE
      END IF

      ! Check for number
      IF (SCAN(c,'0123456789.') > 0) THEN
        r = RealNum (Func(Ind:),ib,in,err)
        IF (err) CALL ParseErrMsg (Ind, FuncStr, 'Invalid number forma&
            &t:  '//Func(Ind+ib-1:Ind+in-2))
        Ind = Ind+in-1
        IF (Ind > lFunc) EXIT
        c = Func(Ind:Ind)
      ELSEIF (c == '_') THEN
        ! Check for constant
        n = ConstantIndex (Func(Ind:))
        IF (n == 0) CALL ParseErrMsg (Ind, FuncStr, 'Invalid constant')
        Ind = Ind+LEN_TRIM(Consts(n))+1
        IF (Ind > lFunc) EXIT
        c=Func(Ind:Ind)
      ELSE
        ! Check for variable
        n = VariableIndex (Func(Ind:),Var,ib,in)
        IF (n == 0) CALL ParseErrMsg (Ind, FuncStr, 'Invalid element: &
            &'//Func(Ind+ib-1:Ind+in-2))
        Ind = Ind+in-1
        IF (Ind > lFunc) EXIT
        c = Func(Ind:Ind)
      END IF

      ! Check for closing parenthesis
      DO WHILE (c == ')')
        IF (.NOT.stack_isempty(functionParenthDepth)) THEN
          IF(stack_backInt(functionParenthDepth) == ParCnt) &
              idummy=stack_popbackInt(functionParenthDepth)
        END IF
        ParCnt = ParCnt-1
        IF (ParCnt < 0) CALL ParseErrMsg (Ind, FuncStr, 'Mismatched pa&
            &renthesis')
        IF (Func(Ind-1:Ind-1) == '(') CALL ParseErrMsg (Ind-1,&
            & FuncStr, 'Empty parentheses')
        Ind = Ind+1
        IF (Ind > lFunc) EXIT
        c = Func(Ind:Ind)
      END DO

      ! Now, we have a legal operand: A legal operator or end of string must follow
999   IF (Ind > lFunc) EXIT
      
      ! Check operators
      opSize = 0
      IF (.NOT.stack_isempty(functionParenthDepth)) THEN
        IF (c == ',' .AND. &
            stack_backInt(functionParenthDepth) == parCnt) THEN
          opSize = 1
        ELSE
          opSize = isOperator(Func(Ind:))
        END IF
      ELSE
        opSize = isOperator(Func(Ind:))
      END IF
      IF (opSize == 0) CALL ParseErrMsg (Ind+1, FuncStr, 'Operator exp&
          &ected')
      
      ! Now, we have an operand and an operator: the next loop will check for another 
      ! operand (must appear)
      Ind = Ind+opSize
    END DO
    IF (ParCnt > 0) CALL ParseErrMsg (Ind, FuncStr, 'Missing )')
    CALL stack_release(functionParenthDepth)
  END SUBROUTINE CheckSyntax

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE ParseErrMsg (j, FuncStr, Msg)

!<description>
    ! Print error message and terminate program
!</description>
    
!<input>
    
    ! Position in function string
    INTEGER, INTENT(in) :: j
    
    ! Original function string
    CHARACTER (LEN=*), INTENT(in) :: FuncStr

    ! OPTIONAL: Error message
    CHARACTER (LEN=*), OPTIONAL, INTENT(in) :: Msg
!</input>
!</subroutine>

    ! local variables
    INTEGER :: k
    
    IF (PRESENT(Msg)) THEN
      WRITE(*,*) '*** Error in syntax of function string: '//Msg
    ELSE
      WRITE(*,*) '*** Error in syntax of function string:'
    ENDIF
    WRITE(*,*)
    WRITE(*,'(A)') ' '//FuncStr
    ! Advance to the jth position
    DO k=1,ipos(j)
      WRITE(*,'(A)',ADVANCE='NO') ' '
    END DO
    WRITE(*,'(A)') '?'
    STOP
  END SUBROUTINE ParseErrMsg

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
      CALL LowCase (str(1:k), fun)
      IF (fun == Funcs(j)) THEN
        ! Found a matching function
        n = j
        EXIT
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
    
    ! Check all math functions
    n = 0
    DO j=C_PI,C_EXP
      k = MIN(LEN_TRIM(Consts(j)), LEN(str(2:)))
      ! Compare lower case letters
      CALL LowCase (str(2:k+1), con)
      IF (con == Consts(j)) THEN
        ! Found a matching constant
        n = j
        EXIT
      END IF
    END DO
  END FUNCTION ConstantIndex

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
    ipos = (/ (k,k=1,lstr) /)
    k = 1
    DO WHILE (str(k:lstr) /= ' ')                             
      IF (str(k:k) == ' ') THEN
        str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
        ipos(k:lstr) = (/ ipos(k+1:lstr), 0 /)             ! Move 1 element to left
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

  SUBROUTINE Compile (Comp, F, Var)

!<description>
    ! Compile i-th function string F into bytecode
!</description>

!<input>
    ! Function string
    CHARACTER (LEN=*), INTENT(in) :: F

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_fparserComponent), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(is), DIMENSION(:), POINTER :: ByteCode
    REAL(DP), DIMENSION(:), POINTER :: Immed
    INTEGER :: ind
    
    IF (ASSOCIATED(Comp%ByteCode)) DEALLOCATE (Comp%ByteCode)
    IF (ASSOCIATED(Comp%Immed))    DEALLOCATE (Comp%Immed)
    Comp%ByteCodeSize = 0
    Comp%ImmedSize    = 0
    Comp%StackSize    = 0
    Comp%StackPtr     = 0

    ! Allocate initial stacks
    ALLOCATE(Comp%ByteCode(1024),Comp%Immed(1024))

    ! Compile string into bytecode
    ind = CompileExpression(Comp,F,1,Var)

    ! Adjust stack sizes
    ALLOCATE(ByteCode(Comp%ByteCodeSize)); ByteCode=Comp%ByteCode
    DEALLOCATE(Comp%ByteCode); ALLOCATE(Comp%ByteCode(Comp%ByteCodeSize))
    Comp%ByteCode=ByteCode; DEALLOCATE(ByteCode)
    
    ALLOCATE(Immed(Comp%ImmedSize)); Immed=Comp%Immed
    DEALLOCATE(Comp%Immed); ALLOCATE(Comp%Immed(Comp%ImmedSize))
    Comp%Immed=Immed; DEALLOCATE(Immed)
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
    IF (Comp%StackPtr > Comp%StackSize)&
        Comp%StackSize=Comp%StackSize+1
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
    INTEGER(is), DIMENSION(:), POINTER :: ByteCode

    Comp%ByteCodeSize = Comp%ByteCodeSize + 1
    IF (ASSOCIATED(Comp%ByteCode)) THEN
      IF (Comp%ByteCodeSize > SIZE(Comp%ByteCode)) THEN
        ALLOCATE(ByteCode(Comp%ByteCodeSize-1))
        ByteCode=Comp%ByteCode
        DEALLOCATE(Comp%ByteCode)
        ALLOCATE(Comp%ByteCode(2*Comp%ByteCodeSize))
        Comp%ByteCode(1:Comp%ByteCodeSize-1)=ByteCode
        DEALLOCATE(ByteCode)
      END IF
      Comp%ByteCode(Comp%ByteCodeSize) = byte
    END IF
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
   
    IF (ASSOCIATED(Comp%ByteCode)) Comp%ByteCode(Comp%ByteCodeSize) = 0
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
   
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: Immed
    
    Comp%ImmedSize = Comp%ImmedSize + 1
    IF (ASSOCIATED(Comp%Immed)) THEN
      IF (Comp%ImmedSize > SIZE(Comp%Immed)) THEN
        ALLOCATE(Immed(Comp%ImmedSize-1))
        Immed=Comp%Immed
        DEALLOCATE(Comp%Immed)
        ALLOCATE(Comp%Immed(2*Comp%ImmedSize))
        Comp%Immed(1:Comp%ImmedSize-1)=Immed
        DEALLOCATE(Immed)
      END IF
      Comp%Immed(Comp%ImmedSize) = immediate

    END IF
  END SUBROUTINE AddImmediate

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
    
    ! local variables
    INTEGER(is), DIMENSION(:), POINTER :: ByteCode

    IF (Comp%useDegreeConversion) THEN
      SELECT CASE(opCode)
      CASE(cCos,Ccosh,cCot,cCsc,cSec,cSin,cSinh,cTan,cTanh)
        CALL AddCompiledByte(Comp,cRad)
      END SELECT
    END IF

    Comp%ByteCodeSize = Comp%ByteCodeSize + 1
    IF (ASSOCIATED(Comp%ByteCode)) THEN
      IF (Comp%ByteCodeSize > SIZE(Comp%ByteCode)) THEN
        ALLOCATE(ByteCode(Comp%ByteCodeSize-1))
        ByteCode=Comp%ByteCode
        DEALLOCATE(Comp%ByteCode)
        ALLOCATE(Comp%ByteCode(2*Comp%ByteCodeSize))
        Comp%ByteCode(1:Comp%ByteCodeSize-1)=ByteCode
        DEALLOCATE(ByteCode)
      END IF
      Comp%ByteCode(Comp%ByteCodeSize) = opcode
    END IF

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
          EXIT ! - otherwise STOP
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
          EXIT ! - otherwise STOP
        END IF
      CASE ('e','E','d','D') ! Permitted only
        IF (InMan) THEN
          Eflag=.true.; InMan=.false. ! - following mantissa
        ELSE
          EXIT ! - otherwise STOP
        ENDIF
      CASE DEFAULT
        EXIT ! STOP at all other characters
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

!<subroutine>

  SUBROUTINE LowCase (str1, str2)

!<description>
    ! Transform upper case letters in str1 into lower case 
    ! letters, result is str2
!</description>
    
!<input>
    ! Source string
    CHARACTER (LEN=*),  INTENT(in) :: str1
!</input>

!<output>
    ! Destination string
    CHARACTER (LEN=*), INTENT(out) :: str2
!</output>
!</subroutine>
    
    ! local variables
    INTEGER :: j,k
    CHARACTER (LEN=*),   PARAMETER :: lc = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    str2 = str1
    DO j=1,LEN_TRIM(str1)
      k = INDEX(uc,str1(j:j))
      IF (k > 0) str2(j:j) = lc(k:k)
    END DO
  END SUBROUTINE LowCase

  ! *****************************************************************************

!<function>

  RECURSIVE FUNCTION CompileExpression(Comp, Func, ind, Var,&
      & stopAtComma) RESULT(ind2)

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
    INTEGER :: lFunc
    
    lFunc=LEN_TRIM(Func)
    
    ind2 = CompileOr(Comp,Func, ind, Var)
    IF(ind2 > lFunc) RETURN

    IF (PRESENT(stopAtComma)) THEN
      IF (stopAtComma) RETURN
    END IF
    
    DO WHILE (Func(ind2:ind2) == ',')
      ind2 = CompileOr(Comp, Func, ind2+1, Var)
      IF (ind2 > lFunc) RETURN
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
    INTEGER :: lFunc

    lFunc=LEN_TRIM(Func)

    ind2 = CompileAnd(Comp, Func, ind, Var)
    IF (ind2 > lFunc) RETURN

    DO WHILE(Func(ind2:ind2) == '|')
      ind2 = CompileAnd(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, cOr)
      Comp%StackPtr = Comp%StackPtr-1
      IF (ind2 > lFunc) RETURN
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
    INTEGER :: lFunc

    lFunc=LEN_TRIM(Func)
    
    ind2 = CompileComparison(Comp, Func, ind, Var)
    IF (ind2 > lFunc) RETURN
    
    DO WHILE(Func(ind2:ind2) == '&')
      ind2 = CompileComparison(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, cAnd)
      Comp%StackPtr = Comp%StackPtr-1
      IF (ind2 > lFunc) RETURN
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
    INTEGER :: lFunc

    lFunc=LEN_TRIM(Func)
    
    ind2 = CompileAddition(Comp, Func, ind, Var)
    IF (ind2 > lFunc) RETURN
    
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
      
      IF (ind2 > lFunc) RETURN
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
    INTEGER :: lFunc

    lFunc=LEN_TRIM(Func)
    
    ind2 = CompileMult(Comp, Func, ind, Var)
    IF (ind2 > lFunc) RETURN
    
    c=Func(ind2:ind2)
    DO WHILE(c == '+' .OR. c == '-')
      ind2 = CompileMult(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, MERGE(cAdd,cSub,c == '+'))
      Comp%StackPtr = Comp%StackPtr-1
      
      IF (ind2 > lFunc) RETURN
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
    INTEGER :: lFunc

    lFunc=LEN_TRIM(Func)
    
    ind2 = CompileUnaryMinus(Comp, Func, ind, Var)
    IF (ind2 > lFunc) RETURN
    
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

      IF (ind2 > lFunc) RETURN
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
    INTEGER :: lFunc

    lFunc=LEN_TRIM(Func)

    c = Func(ind:ind)
    IF (c == '-' .OR. c == '!') THEN
      ind2 = ind+1
      IF (ind2 > lFunc) RETURN
      ind2 = CompilePow(Comp, Func, ind2, Var)

      ! If we are negating a constant, negate the constant itself
      IF (c == '-' .AND. &
          & Comp%ByteCode(Comp%ByteCodeSize) == cImmed) THEN
        Comp%Immed(Comp%ImmedSize) = -Comp%Immed(Comp%ImmedSize)
        
        ! If we are negating a negation, we can remove both
      ELSEIF (c == '-' .AND. &
          & Comp%ByteCode(Comp%ByteCodeSize) == cNeg) THEN
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
    INTEGER :: lFunc

    lFunc=LEN_TRIM(Func)

    ind2 = CompileElement(Comp, Func, ind, Var)

    IF (ind2 > lFunc) RETURN

    DO WHILE(Func(ind2:ind2) == '^')
      ind2 = CompileUnaryMinus(Comp, Func, ind2+1, Var)

      CALL AddCompiledByte(Comp, cPow)
      Comp%StackPtr = Comp%StackPtr-1
      IF (ind2 > lFunc) RETURN
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
    REAL(DP) :: r
    INTEGER(is) :: n
    INTEGER :: ind1,ib,in,requiredParams
    LOGICAL :: err

    ind1=ind; c=Func(ind1:ind1)
    IF (c == '(') THEN
      ind1 = CompileExpression(Comp, Func, ind1+1, Var, .FALSE.)
      ind2 = ind1+1   ! Func(ind1:ind1) is ')'
      RETURN
    END IF

    ! Check for numbers
    IF (SCAN(c,'0123456789,') > 0) THEN
      r = RealNum (Func(ind1:),ib,in,err)
      IF (err) CALL ParseErrMsg (ind1, Func, 'Invalid number format: &
          & '//Func(ind1+ib-1:ind1+in-2))
      CALL AddImmediate(Comp, r)
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
          RETURN        
        END IF

        requiredParams=MathFunctionParameters(n)        
        ind2 = CompileFunctionParameters(Comp, Func, ind2+1, Var, &
            & requiredParams)
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

      ! Check for variables
      n = VariableIndex(Func(ind1:), Var, ib, in)
      IF (n > 0) n = VarBegin+n-1
      CALL AddCompiledByte(Comp, n)
      CALL incStackPtr(Comp)
      ind2 = ind1+in-1
      RETURN
    END IF

    CALL ParseErrMsg (ind1, Func, 'An unexpected error occured. Please&
        & make a full bug report to the author')
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

      IF (Comp%StackPtr /= curStackPtr+requiredParams)&
          &CALL ParseErrMsg (ind, Func, 'Illegal number of parameters &
          &to function')

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
    INTEGER :: lFunc,curByteCodeSize,curByteCodeSize2,curImmedSize2

    lFunc=LEN_TRIM(Func)

    ind2 = CompileExpression(Comp, Func, ind, Var, .TRUE.) ! Condition branch
    IF (ind2 > lFunc) RETURN

    IF (Func(ind2:ind2) /= ',') CALL ParseErrMsg (ind2, Func, 'Illegal&
        & number of parameters to function')
    CALL AddCompiledByte(Comp, cIf)
    curByteCodeSize = Comp%ByteCodeSize
    CALL AddCompiledByte(Comp, 0_is) ! Jump index will be set below
    CALL AddCompiledByte(Comp, 0_is) ! Immed jump index will be set below
    Comp%StackPtr = Comp%StackPtr-1

    ind2 = CompileExpression(Comp, Func, ind2+1, Var, .TRUE.) ! Then branch
    IF (ind2 > lFunc) RETURN

    IF (Func(ind2:ind2) /= ',') CALL ParseErrMsg (ind2, Func, 'Illegal&
        & number of parameters to function')
    CALL AddCompiledByte(Comp, cJump)
    curByteCodeSize2 = Comp%ByteCodeSize
    curImmedSize2 = Comp%ImmedSize
    CALL AddCompiledByte(Comp, 0_is) ! Jump index will be set below
    CALL AddCompiledByte(Comp, 0_is) ! Immed jump index will be set below
    Comp%StackPtr = Comp%StackPtr-1

    ind2 = CompileExpression(Comp, Func, ind2+1, Var, .TRUE.) ! Else branch
    IF (ind2 > lFunc) RETURN

    IF (Func(ind2:ind2) /= ')') CALL ParseErrMsg (ind2, Func, 'Illegal&
        & number of parameters to function')
    
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
      
    r=d*(Constvals(C_PI)/180._DP)
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
    
    d=r*(180._DP/Constvals(c_PI))
  END FUNCTION RadToDeg
END MODULE fparser
