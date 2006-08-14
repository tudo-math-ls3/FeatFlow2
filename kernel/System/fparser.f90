!##############################################################################
!# ****************************************************************************
!# <name> fparser </name>
!# ****************************************************************************
!#
!# <purpose>
!# This public domain function parser module is intended for applications 
!# where a set of mathematical expressions is specified at runtime and is 
!# then evaluated for a large number of variable values. This is done by
!# compiling the set of function strings into byte code, which is interpreted
!# very efficiently for the various variable values. 
!#
!# The source code is available from:
!# http://www.its.uni-karlsruhe.de/~schmehl/opensource/fparser-v1.1.tar.gz
!#
!# Please send comments, corrections or questions to the author:
!# Roland Schmehl <Roland.Schmehl@mach.uni-karlsruhe.de>
!#
!#
!# The function parser concept is based on a C++ class library written by Warp 
!# <warp@iki.fi> available from:
!# http://www.students.tut.fi/~warp/FunctionParser/fparser.zip
!#
!# 
!# The following description of this module is taken from the
!# original README file
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
!# This command imports only 5 public names: initf, parsef, evalf, EvalErrMsg
!# and EvalErrType, which are explained in the following. The remainder of the 
!# module is hidden to the calling program.
!#
!# Step 1 - Initialization
!# -----------------------
!# The parser module has to be initialized for the simultaneous evaluation of 
!# n functions by calling the module subroutine initp one time in your Fortran 
!# code:
!# 
!# CALL initf (n)
!# 
!# This allocates i=1,...,n internal data structures used by the byte-compiler 
!# and subsequently by the bytecode-interpreter.
!#
!# Step 2 - Function parsing
!# -------------------------
!# The i-th function string FuncStr is parsed (checked and compiled) into the 
!# i-th bytecode by calling the module subroutine parsef:
!#
!# CALL fparser_parseFunction (i, FuncStr, Var)
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
!# a = fparser_evalFunction (i, Val)
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
!# sets the error flag EvalErrType (also imported by the USE statement) to 
!# a value > 0 (EvalErrType=0 indicates no error). An error message from the 
!# bytecode-interpreter can be obtained by calling the character function 
!# fparser_ErrorMsg () without any argument.
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
!# Mathematical Operators recognized are +, -, *, /, ** or alternatively ^, 
!# whereas symbols for brackets must be (). 
!# In addition the the min and max functions for two(!) arguments
!# have been implemented
!#
!# The function parser recognizes the (single argument) Fortran 90 intrinsic 
!# functions abs, exp, log10, log, sqrt, sinh, cosh, tanh, sin, cos, tan, asin, 
!# acos, atan. Parsing for intrinsic functions is case INsensitive.
!#
!# Operations are evaluated in the correct order:
!#
!#  ()          expressions in brackets first
!#  -A          unary minus (or plus)
!#  A**B A^B    exponentiation (A raised to the power B)
!#  A*B  A/B    multiplication and division
!#  A+B  A-B    addition and subtraction
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
!# 4.) fparser_evalFunction
!#     -> Evaluate bytecode
!#
!# 5.) fparser_ErrorMsg
!#     -> Get error message from function parser
!#
!# </purpose>
!##############################################################################

MODULE fparser
  USE fsystem
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_Comp
  PUBLIC :: fparser_create
  PUBLIC :: fparser_release
  PUBLIC :: fparser_parseFunction
  PUBLIC :: fparser_evalFunction
  PUBLIC :: fparser_ErrorMsg

!<constants>

!<constantblock description="data type for parser bytecode">
  ! Data type of bytecode
  INTEGER, PARAMETER :: is = SELECTED_INT_KIND(1)
!</constantblock>

!<constantblock description="keywords for parser">
  INTEGER(is), PARAMETER :: cImmed   = 1, &
                            cNeg     = 2, &
                            cAdd     = 3, & 
                            cSub     = 4, & 
                            cMul     = 5, & 
                            cDiv     = 6, & 
                            cPow     = 7, & 
                            cMin     = 8, &
                            cMax     = 9, &
                            cAbs     = 10,&
                            cExp     = 11,&
                            cLog10   = 12,&
                            cLog     = 13,&
                            cSqrt    = 14,&
                            cSinh    = 15,&
                            cCosh    = 16,&
                            cTanh    = 17,&
                            cSin     = 18,&
                            cCos     = 19,&
                            cTan     = 20,&
                            cAsin    = 21,&
                            cAcos    = 22,&
                            cAtan    = 23,&
                            VarBegin = 24
!</constantblock>

!<constantblock description="operands for parser">
  CHARACTER (LEN=1), DIMENSION(cAdd:cMax),  PARAMETER :: Ops      = (/ '+',     &
                                                                       '-',     &
                                                                       '*',     &
                                                                       '/',     &
                                                                       '^',     &
                                                                       'm',     &
                                                                       'M' /)
!</constantblock>

!<constantblock description="functionss for parser">
  CHARACTER (LEN=5), DIMENSION(cAbs:cAtan), PARAMETER :: Funcs    = (/ 'abs  ', &
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
                                                                       'asin ', &
                                                                       'acos ', &
                                                                       'atan ' /)
!</constantblock>

!</constants>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<globals>

  ! Evaluation error code: =0: no error occured, >0
  INTEGER, PUBLIC :: EvalErrType 

  ! Associates function strings
  INTEGER,   DIMENSION(:),  ALLOCATABLE :: ipos
!</globals>

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<types>

!<typeblock>

  ! Type block for storing the bytecode of the function parser
  TYPE t_Comp
    ! Bytecode
    INTEGER(is), DIMENSION(:), POINTER :: ByteCode
    
    ! Size of bytecode
    INTEGER :: ByteCodeSize
    
    !
    REAL(DP), DIMENSION(:), POINTER :: Immed
    
    !
    INTEGER :: ImmedSize
    
    ! Stack
    REAL(DP), DIMENSION(:), POINTER :: Stack
    
    ! Stack size
    INTEGER :: StackSize
    
    ! Stack pointer
    INTEGER :: StackPtr
  END TYPE t_Comp
!</typeblock>


!</types>


CONTAINS
  
  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE fparser_create (Comp,n)

!<description>
    ! Initialize function parser for n functions
!</description>

!<input>
    ! Number of functions
    INTEGER, INTENT(in) :: n                                 
!</input>

!<output>
    ! Function parser object
    TYPE(t_Comp), DIMENSION(:), POINTER :: Comp
!</output>
!</subroutine>
    
    ! local variables
    INTEGER :: i

    ALLOCATE (Comp(n))
    DO i=1,n
      NULLIFY (Comp(i)%ByteCode,Comp(i)%Immed,Comp(i)%Stack)
    END DO
  END SUBROUTINE fparser_create

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE fparser_release (Comp)

!<description>
    ! Release function parser
!</description>

!<input>
    ! Function parser
    TYPE(t_Comp), DIMENSION(:), POINTER :: Comp
!</input>
!</subroutine>
    
    ! local variables
    INTEGER :: i
    
    DO i=1,SIZE(Comp)
      DEALLOCATE(Comp(i)%ByteCode,Comp(i)%Immed,Comp(i)%Stack)
    END DO
    DEALLOCATE(Comp)
  END SUBROUTINE fparser_release
  
  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE fparser_parseFunction (Comp, i, FuncStr, Var)

!<description>
    ! Parse ith function string FuncStr and compile it into bytecode
!</description>
    
!<input>
    ! Function identifier
    INTEGER, INTENT(in) :: i

    ! Function string
    CHARACTER (LEN=*), INTENT(in) :: FuncStr

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_Comp),  DIMENSION(:), INTENT(inout) :: Comp
!</inputoutput>
!</subroutine>

    ! local variables
    CHARACTER (LEN=LEN(FuncStr)) :: Func
    
    IF (i < 1 .OR. i > SIZE(Comp)) THEN
      WRITE(*,*) '*** Parser error: Function number ',i,' out of range'
      STOP
    END IF
    ALLOCATE (ipos(LEN_TRIM(FuncStr)))                       ! Char. positions in orig. string
    Func = FuncStr                                           ! Local copy of function string
    CALL Replace ('**','^ ',Func)                            ! Exponent into 1-Char. format
    CALL Replace ('MIN','m  ',Func)
    CALL Replace ('min','m  ',Func)
    CALL Replace ('MAX','M  ',Func)
    CALL Replace ('max','M  ',Func)
    CALL RemoveSpaces (Func)                                 ! Condense function string
    CALL CheckSyntax (Func,FuncStr,Var)
    DEALLOCATE (ipos)
    CALL Compile (Comp, i,Func,Var)                                ! Compile into bytecode
  END SUBROUTINE fparser_parseFunction

  ! *****************************************************************************

!<function>

  FUNCTION fparser_evalFunction (Comp, i, Val) RESULT (res)

!<description>
    ! Evaluate bytecode of ith function for the values passed in
    ! array Val(:)
!</description>

!<input>
    ! Function parser
    TYPE (t_Comp),  DIMENSION(:), INTENT(in) :: Comp

    ! Function identifier
    INTEGER, INTENT(in) :: i

    ! Variable values
    REAL(DP), DIMENSION(:), INTENT(in) :: Val
!</input>

!<result>
    ! Evaluated function
    REAL(DP) :: res
!</result>
!</function>
    
    ! local variables
    INTEGER :: IP, DataP, StackP
    REAL(DP), PARAMETER :: zero = 0._DP
    
    DataP  = 1
    StackP = 0
    DO IP=1,Comp(i)%ByteCodeSize
      SELECT CASE (Comp(i)%ByteCode(IP))
        
      CASE (cImmed); StackP=StackP+1; Comp(i)%Stack(StackP)=Comp(i)%Immed(DataP); DataP=DataP+1
      CASE   (cNeg); Comp(i)%Stack(StackP)=-Comp(i)%Stack(StackP)
      CASE   (cAdd); Comp(i)%Stack(StackP-1)=Comp(i)%Stack(StackP-1)+Comp(i)%Stack(StackP); StackP=StackP-1
      CASE   (cSub); Comp(i)%Stack(StackP-1)=Comp(i)%Stack(StackP-1)-Comp(i)%Stack(StackP); StackP=StackP-1
      CASE   (cMul); Comp(i)%Stack(StackP-1)=Comp(i)%Stack(StackP-1)*Comp(i)%Stack(StackP); StackP=StackP-1
      CASE   (cDiv); IF (Comp(i)%Stack(StackP)==0._DP) THEN; EvalErrType=1; res=zero; RETURN; ENDIF
        Comp(i)%Stack(StackP-1)=Comp(i)%Stack(StackP-1)/Comp(i)%Stack(StackP); StackP=StackP-1
      CASE   (cPow); Comp(i)%Stack(StackP-1)=Comp(i)%Stack(StackP-1)**Comp(i)%Stack(StackP); StackP=StackP-1
      CASE   (cMin); Comp(i)%Stack(StackP-1)=MIN(Comp(i)%Stack(StackP-1),Comp(i)%Stack(StackP)); StackP=StackP-1
      CASE   (cMax); Comp(i)%Stack(StackP-1)=MAX(Comp(i)%Stack(StackP-1),Comp(i)%Stack(StackP)); StackP=StackP-1
      CASE   (cAbs); Comp(i)%Stack(StackP)=ABS(Comp(i)%Stack(StackP))
      CASE   (cExp); Comp(i)%Stack(StackP)=EXP(Comp(i)%Stack(StackP))
      CASE (cLog10); IF (Comp(i)%Stack(StackP)<=0._DP) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
        Comp(i)%Stack(StackP)=LOG10(Comp(i)%Stack(StackP))
      CASE   (cLog); IF (Comp(i)%Stack(StackP)<=0._DP) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
        Comp(i)%Stack(StackP)=LOG(Comp(i)%Stack(StackP))
      CASE  (cSqrt); IF (Comp(i)%Stack(StackP)<0._DP) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
        Comp(i)%Stack(StackP)=SQRT(Comp(i)%Stack(StackP))
      CASE  (cSinh); Comp(i)%Stack(StackP)=SINH(Comp(i)%Stack(StackP))
      CASE  (cCosh); Comp(i)%Stack(StackP)=COSH(Comp(i)%Stack(StackP))
      CASE  (cTanh); Comp(i)%Stack(StackP)=TANH(Comp(i)%Stack(StackP))
      CASE   (cSin); Comp(i)%Stack(StackP)=SIN(Comp(i)%Stack(StackP))
      CASE   (cCos); Comp(i)%Stack(StackP)=COS(Comp(i)%Stack(StackP))
      CASE   (cTan); Comp(i)%Stack(StackP)=TAN(Comp(i)%Stack(StackP))
      CASE  (cAsin); IF ((Comp(i)%Stack(StackP)<-1._DP).OR.(Comp(i)%Stack(StackP)>1._DP)) THEN
        EvalErrType=4; res=zero; RETURN; ENDIF
        Comp(i)%Stack(StackP)=ASIN(Comp(i)%Stack(StackP))
      CASE  (cAcos); IF ((Comp(i)%Stack(StackP)<-1._DP).OR.(Comp(i)%Stack(StackP)>1._DP)) THEN
        EvalErrType=4; res=zero; RETURN; ENDIF
        Comp(i)%Stack(StackP)=ACOS(Comp(i)%Stack(StackP))
      CASE  (cAtan); Comp(i)%Stack(StackP)=ATAN(Comp(i)%Stack(StackP))
      CASE  DEFAULT; StackP=StackP+1; Comp(i)%Stack(StackP)=Val(Comp(i)%ByteCode(IP)-VarBegin+1)
      END SELECT
    END DO
    EvalErrType = 0
    res = Comp(i)%Stack(1)
  END FUNCTION fparser_evalFunction

  ! *****************************************************************************

!<function>

  FUNCTION fparser_ErrorMsg () RESULT (msg)

!<description>
    ! Return error message of function parser
!</description>

    ! local constants
    CHARACTER (LEN=*), DIMENSION(4), PARAMETER :: m = (/ 'Division by zero                ', &
                                                         'Argument of SQRT negative       ', &
                                                         'Argument of LOG negative        ', &
                                                         'Argument of ASIN or ACOS illegal' /)
    
!<result>
    ! Error messages
    CHARACTER (LEN=LEN(m)) :: msg
!>/result>
!</function>

    
    
    IF (EvalErrType < 1 .OR. EvalErrType > SIZE(m)) THEN
      msg = ''
    ELSE
      msg = m(EvalErrType)
    ENDIF
  END FUNCTION fparser_ErrorMsg

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE CheckSyntax (Func,FuncStr,Var)

!<description>
    ! Check syntax of function string,  returns 0 if syntax is ok
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
    INTEGER(is) :: n
    CHARACTER (LEN=1) :: c
    REAL(DP) :: r
    LOGICAL :: err
    INTEGER :: ParCnt,j,ib,in,lFunc
    
    j = 1
    ParCnt = 0
    lFunc = LEN_TRIM(Func)
    step: DO
      IF (j > lFunc) CALL ParseErrMsg (j, FuncStr)
      c = Func(j:j)
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      ! Check for valid operand (must appear)
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      IF (c == '-' .OR. c == '+') THEN                      ! Check for leading - or +
        j = j+1
        IF (j > lFunc) CALL ParseErrMsg (j, FuncStr, 'Missing operand')
        c = Func(j:j)
        IF (ANY(c == Ops)) CALL ParseErrMsg (j, FuncStr, 'Multiple operators')
      END IF
      n = MathFunctionIndex (Func(j:))
      IF (n > 0) THEN                                       ! Check for math function
        j = j+LEN_TRIM(Funcs(n))
        IF (j > lFunc) CALL ParseErrMsg (j, FuncStr, 'Missing function argument')
        c = Func(j:j)
        IF (c /= '(') CALL ParseErrMsg (j, FuncStr, 'Missing opening parenthesis')
      END IF
      IF (c == '(') THEN                                    ! Check for opening parenthesis
        ParCnt = ParCnt+1
        j = j+1
        CYCLE step
      END IF
      IF (SCAN(c,'0123456789.') > 0) THEN                   ! Check for number
        r = RealNum (Func(j:),ib,in,err)
        IF (err) CALL ParseErrMsg (j, FuncStr, 'Invalid number format:  '//Func(j+ib-1:j+in-2))
        j = j+in-1
        IF (j > lFunc) EXIT
        c = Func(j:j)
      ELSE                                                  ! Check for variable
        n = VariableIndex (Func(j:),Var,ib,in)
        IF (n == 0) CALL ParseErrMsg (j, FuncStr, 'Invalid element: '//Func(j+ib-1:j+in-2))
        j = j+in-1
        IF (j > lFunc) EXIT
        c = Func(j:j)
      END IF
      DO WHILE (c == ')')                                   ! Check for closing parenthesis
        ParCnt = ParCnt-1
        IF (ParCnt < 0) CALL ParseErrMsg (j, FuncStr, 'Mismatched parenthesis')
        IF (Func(j-1:j-1) == '(') CALL ParseErrMsg (j-1, FuncStr, 'Empty parentheses')
        j = j+1
        IF (j > lFunc) EXIT
        c = Func(j:j)
      END DO
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      ! Now, we have a legal operand: A legal operator or end of string must follow
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      IF (j > lFunc) EXIT
      IF (ANY(c == Ops)) THEN                               ! Check for multiple operators
        IF (j+1 > lFunc) CALL ParseErrMsg (j, FuncStr)
        IF (ANY(Func(j+1:j+1) == Ops)) CALL ParseErrMsg (j+1, FuncStr, 'Multiple operators')
      ELSE                                                  ! Check for next operand
        CALL ParseErrMsg (j, FuncStr, 'Missing operator')
      END IF
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      ! Now, we have an operand and an operator: the next loop will check for another 
      ! operand (must appear)
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      j = j+1
    END DO step
    IF (ParCnt > 0) CALL ParseErrMsg (j, FuncStr, 'Missing )')
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
    DO k=1,ipos(j)
      WRITE(*,'(A)',ADVANCE='NO') ' '                       ! Advance to the jth position
    END DO
    WRITE(*,'(A)') '?'
    STOP
  END SUBROUTINE ParseErrMsg

  ! *****************************************************************************

!<function>

  FUNCTION OperatorIndex (c) RESULT (n)

!<description>
    ! Return operator index
!</description>

!<input>
    ! Operator string
    CHARACTER (LEN=1), INTENT(in) :: c
!</input>

!<result>
    INTEGER(is) :: n
!</result>
!</function>
    
    ! local variables
    INTEGER(iS) :: j
    
    n = 0
    DO j=cAdd,cMax
      IF (c == Ops(j)) THEN
        n = j
        EXIT
      END IF
    END DO
  END FUNCTION OperatorIndex

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
    
    n = 0
    DO j=cAbs,cAtan                                          ! Check all math functions
      k = MIN(LEN_TRIM(Funcs(j)), LEN(str))   
      CALL LowCase (str(1:k), fun)
      IF (fun == Funcs(j)) THEN                             ! Compare lower case letters
        n = j                                              ! Found a matching function
        EXIT
      END IF
    END DO
  END FUNCTION MathFunctionIndex

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
      DO ib=1,lstr                                          ! Search for first character in str
        IF (str(ib:ib) /= ' ') EXIT                        ! When lstr>0 at least 1 char in str
      END DO
      DO in=ib,lstr                                         ! Search for name terminators
        IF (SCAN(str(in:in),'+-*/^m) ') > 0) EXIT
      END DO
      DO j=1,SIZE(Var)
        IF (str(ib:in-1) == Var(j)) THEN                     
          n = j                                           ! Variable name found
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

  SUBROUTINE Compile (Comp, i, F, Var)

!<description>
    ! Compile i-th function string F into bytecode
!</description>

!<input>
    ! Function identifier
    INTEGER, INTENT(in) :: i

    ! Function string
    CHARACTER (LEN=*), INTENT(in) :: F

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_Comp), DIMENSION(:), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: istat
    
    IF (ASSOCIATED(Comp(i)%ByteCode)) DEALLOCATE ( Comp(i)%ByteCode, &
                                                   Comp(i)%Immed,    &
                                                   Comp(i)%Stack     )
    Comp(i)%ByteCodeSize = 0
    Comp(i)%ImmedSize    = 0
    Comp(i)%StackSize    = 0
    Comp(i)%StackPtr     = 0
    CALL CompileSubstr (Comp,i,F,1,LEN_TRIM(F),Var)               ! Compile string to determine size
    ALLOCATE ( Comp(i)%ByteCode(Comp(i)%ByteCodeSize), & 
               Comp(i)%Immed(Comp(i)%ImmedSize),       &
               Comp(i)%Stack(Comp(i)%StackSize),       &
               STAT = istat                            )
    IF (istat /= 0) THEN
      WRITE(*,*) '*** Parser error: Memmory allocation for byte code failed'
      STOP
    ELSE
      Comp(i)%ByteCodeSize = 0
      Comp(i)%ImmedSize    = 0
      Comp(i)%StackSize    = 0
      Comp(i)%StackPtr     = 0
      CALL CompileSubstr (Comp,i,F,1,LEN_TRIM(F),Var)            ! Compile string into bytecode
    END IF
    !
  END SUBROUTINE Compile

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE AddCompiledByte (Comp, i, b)

!<description>
    ! Add compiled byte to bytecode
!</description>

!<input>
    ! Function identifier  
    INTEGER, INTENT(in) :: i

    ! Value of byte to be added
    INTEGER(is), INTENT(in) :: b
!</input>

!<inputoutput>
    TYPE (t_Comp), DIMENSION(:), INTENT(INOUT) :: Comp
!</inputoutput>
!</subroutine>
   
    Comp(i)%ByteCodeSize = Comp(i)%ByteCodeSize + 1
    IF (ASSOCIATED(Comp(i)%ByteCode)) Comp(i)%ByteCode(Comp(i)%ByteCodeSize) = b
  END SUBROUTINE AddCompiledByte

  ! *****************************************************************************

!<function>

  FUNCTION MathItemIndex (Comp, i, F, Var) RESULT (n)

!<description>
    ! Return math item index, if item is real number, enter it into Comp-structure
!</description>

!<input>
    ! Function identifier  
    INTEGER, INTENT(in) :: i

    ! Function substring
    CHARACTER (LEN=*), INTENT(in) :: F

    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_Comp), DIMENSION(:), INTENT(INOUT) :: Comp
!</inputoutput>

!<result>
    ! Byte value of math item
    INTEGER(is) :: n
!</result>
!</function>
    
    n = 0
    IF (SCAN(F(1:1),'0123456789.') > 0) THEN                 ! Check for begin of a number
      Comp(i)%ImmedSize = Comp(i)%ImmedSize + 1
      IF (ASSOCIATED(Comp(i)%Immed)) Comp(i)%Immed(Comp(i)%ImmedSize) = RealNum (F)
      n = cImmed
    ELSE                                                     ! Check for a variable
      n = VariableIndex (F, Var)
      IF (n > 0) n = VarBegin+n-1
    END IF
  END FUNCTION MathItemIndex

  ! *****************************************************************************

!<function>

  FUNCTION CompletelyEnclosed (F, b, e) RESULT (res)

!<description>
    ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
!</dscription>

!<input>
    ! Function substring
    CHARACTER (LEN=*), INTENT(in) :: F

    ! First position of substring
    INTEGER, INTENT(in) :: b

    ! Last position of substring
    INTEGER, INTENT(in) :: e
!</input>

!<result>
    ! TRUE if function substring F(b:e) is completely enclosed by a
    ! pair of paranthesis
    LOGICAL :: res
!</result>
!</function>
    
    ! local variables
    INTEGER :: j,k
    
    res=.false.
    IF (F(b:b) == '(' .AND. F(e:e) == ')') THEN
      k = 0
      DO j=b+1,e-1
        IF     (F(j:j) == '(') THEN
          k = k+1
        ELSEIF (F(j:j) == ')') THEN
          k = k-1
        END IF
        IF (k < 0) EXIT
      END DO
      IF (k == 0) res=.true.                                ! All opened parenthesis closed
    END IF
  END FUNCTION CompletelyEnclosed

  ! *****************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE CompileSubstr (Comp, i, F, b, e, Var)

!<description>
    ! Compile i-th function string F into bytecode
!</decsription>

!<input>
    ! Function identifier 
    INTEGER, INTENT(in) :: i

    ! Function substring
    CHARACTER (LEN=*), INTENT(in) :: F

    ! Begin position substring
    INTEGER, INTENT(in) :: b

    ! End position substring
    INTEGER, INTENT(in) :: e
    
    ! Array with variable names
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var
!</input>

!<inputoutput>
    ! Function parser
    TYPE (t_Comp),  DIMENSION(:),  INTENT(inout) :: Comp
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(is) :: n
    INTEGER :: b2,j,k,io
    CHARACTER (LEN=*), PARAMETER :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
                                             'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for special cases of substring
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (F(b:b) == '+') THEN                              ! Case 1: F(b:e) = '+...'
      CALL CompileSubstr (Comp, i, F, b+1, e, Var)
      RETURN
    ELSEIF (CompletelyEnclosed (F, b, e)) THEN           ! Case 2: F(b:e) = '(...)'
      CALL CompileSubstr (Comp, i, F, b+1, e-1, Var)
      RETURN
    ELSEIF (SCAN(F(b:b),calpha) > 0) THEN        
      n = MathFunctionIndex (F(b:e))
      IF (n > 0) THEN
        b2 = b+INDEX(F(b:e),'(')-1
        IF (CompletelyEnclosed(F, b2, e)) THEN           ! Case 3: F(b:e) = 'fcn(...)'
          CALL CompileSubstr(Comp, i, F, b2+1, e-1, Var)
          CALL AddCompiledByte (Comp, i, n)
          RETURN
        END IF
      END IF
    ELSEIF (F(b:b) == '-') THEN
      IF (CompletelyEnclosed (F, b+1, e)) THEN           ! Case 4: F(b:e) = '-(...)'
        CALL CompileSubstr (Comp, i, F, b+2, e-1, Var)
        CALL AddCompiledByte (Comp, i, cNeg)
        RETURN
      ELSEIF (SCAN(F(b+1:b+1),calpha) > 0) THEN
        n = MathFunctionIndex (F(b+1:e))
        IF (n > 0) THEN
          b2 = b+INDEX(F(b+1:e),'(')
          IF (CompletelyEnclosed(F, b2, e)) THEN         ! Case 5: F(b:e) = '-fcn(...)'
            CALL CompileSubstr(Comp, i, F, b2+1, e-1, Var)
            CALL AddCompiledByte (Comp, i, n)
            CALL AddCompiledByte (Comp, i, cNeg)
            RETURN
          END IF
        END IF
      ENDIF
    END IF

    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
    !----- -------- --------- --------- --------- --------- --------- --------- -------    
    DO io=cAdd,cMax                                      ! Increasing priority +-*/^
      k = 0
      DO j=e,b,-1
        IF     (F(j:j) == ')') THEN
          k = k+1
        ELSEIF (F(j:j) == '(') THEN
          k = k-1
        END IF
        IF (k == 0 .AND. F(j:j) == Ops(io) .AND. IsBinaryOp (j, F)) THEN
          IF (ANY(F(j:j) == Ops(cMul:cMax)) .AND. F(b:b) == '-') THEN ! Case 6: F(b:e) = '-...Op...' with Op > -
            CALL CompileSubstr (Comp, i, F, b+1, e, Var)
            CALL AddCompiledByte (Comp, i, cNeg)
            RETURN                 
          ELSE                                                        ! Case 7: F(b:e) = '...BinOp...'
            CALL CompileSubstr (Comp, i, F, b, j-1, Var)
            CALL CompileSubstr (Comp, i, F, j+1, e, Var)
            CALL AddCompiledByte (Comp, i, OperatorIndex(Ops(io)))
            Comp(i)%StackPtr = Comp(i)%StackPtr - 1
            RETURN
          END IF
        END IF
      END DO
    END DO

    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for remaining items, i.e. variables or explicit numbers
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    b2 = b
    IF (F(b:b) == '-') b2 = b2+1
    n = MathItemIndex(Comp, i, F(b2:e), Var)
    CALL AddCompiledByte (Comp, i, n)
    Comp(i)%StackPtr = Comp(i)%StackPtr + 1
    IF (Comp(i)%StackPtr > Comp(i)%StackSize) Comp(i)%StackSize = Comp(i)%StackSize + 1
    IF (b2 > b) CALL AddCompiledByte (Comp, i, cNeg)
  END SUBROUTINE CompileSubstr

  ! *****************************************************************************

!<function>

  FUNCTION IsBinaryOp (j, F) RESULT (res)

!<description>
    ! Check if operator F(j:j) in string F is binary operator
    ! Special cases already covered elsewhere:              (that is corrected in v1.1)
    ! - operator character F(j:j) is first character of string (j=1)
!</description>

!<input>
    ! Position of Operator
    INTEGER, INTENT(in) :: j

    ! String
    CHARACTER (LEN=*), INTENT(in) :: F
!</input>

!<result>
    ! TRUE if F(j:j) is binary operator
    LOGICAL :: res
!</result>
!</function>

    ! local variables
    INTEGER :: k
    LOGICAL :: Dflag,Pflag

    res=.true.
    IF (F(j:j) == '+' .OR. F(j:j) == '-') THEN               ! Plus or minus sign:
      IF (j == 1) THEN                                      ! - leading unary operator ?
        res = .false.
      ELSEIF (SCAN(F(j-1:j-1),'+-*/^m(') > 0) THEN           ! - other unary operator ?
        res = .false.
      ELSEIF (SCAN(F(j+1:j+1),'0123456789') > 0 .AND. &     ! - in exponent of real number ?
          SCAN(F(j-1:j-1),'eEdD')       > 0) THEN
        Dflag=.false.; Pflag=.false.
        k = j-1
        DO WHILE (k > 1)                                   !   step to the left in mantissa 
          k = k-1
          IF     (SCAN(F(k:k),'0123456789') > 0) THEN
            Dflag=.true.
          ELSEIF (F(k:k) == '.') THEN
            IF (Pflag) THEN
              EXIT                                      !   * EXIT: 2nd appearance of '.'
            ELSE
              Pflag=.true.                              !   * mark 1st appearance of '.'
            ENDIF
          ELSE
            EXIT                                         !   * all other characters
          END IF
        END DO
        IF (Dflag .AND. (k == 1 .OR. SCAN(F(k:k),'+-*/^m(') > 0)) res = .false.
      END IF
    END IF
  END FUNCTION IsBinaryOp

  ! *****************************************************************************

!<function>

  FUNCTION RealNum (str, ibegin, inext, error) RESULT (res)

!<description>
    ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
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
      CASE (' ')                                            ! Only leading blanks permitted
        ib = ib+1
        IF (InMan .OR. Eflag .OR. InExp) EXIT
      CASE ('+','-')                                        ! Permitted only
        IF     (Bflag) THEN           
          InMan=.true.; Bflag=.false.                     ! - at beginning of mantissa
        ELSEIF (Eflag) THEN               
          InExp=.true.; Eflag=.false.                     ! - at beginning of exponent
        ELSE
          EXIT                                            ! - otherwise STOP
        ENDIF
      CASE ('0':'9')                                        ! Mark
        IF     (Bflag) THEN           
          InMan=.true.; Bflag=.false.                     ! - beginning of mantissa
        ELSEIF (Eflag) THEN               
          InExp=.true.; Eflag=.false.                     ! - beginning of exponent
        ENDIF
        IF (InMan) DInMan=.true.                           ! Mantissa contains digit
        IF (InExp) DInExp=.true.                           ! Exponent contains digit
      CASE ('.')
        IF     (Bflag) THEN
          Pflag=.true.                                    ! - mark 1st appearance of '.'
          InMan=.true.; Bflag=.false.                     !   mark beginning of mantissa
        ELSEIF (InMan .AND..NOT.Pflag) THEN
          Pflag=.true.                                    ! - mark 1st appearance of '.'
        ELSE
          EXIT                                            ! - otherwise STOP
        END IF
      CASE ('e','E','d','D')                                ! Permitted only
        IF (InMan) THEN
          Eflag=.true.; InMan=.false.                     ! - following mantissa
        ELSE
          EXIT                                            ! - otherwise STOP
        ENDIF
      CASE DEFAULT
        EXIT                                               ! STOP at all other characters
      END SELECT
      in = in+1
    END DO
    err = (ib > in-1) .OR. (.NOT.DInMan) .OR. ((Eflag.OR.InExp).AND..NOT.DInExp)
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
    ! Transform upper case letters in str1 into lower case letters, result is str2
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
END MODULE fparser
