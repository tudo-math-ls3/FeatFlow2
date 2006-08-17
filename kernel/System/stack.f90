!##############################################################################
!# ****************************************************************************
!# <name> stack </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a dynamic stack which can be use to store
!# integer, single and double data.
!#
!# </purpose>
!##############################################################################

MODULE stack
  USE storage
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_stack
  PUBLIC :: stack_create
  PUBLIC :: stack_release
  PUBLIC :: stack_clear
  PUBLIC :: stack_isempty
  PUBLIC :: stack_size
  PUBLIC :: stack_pushback
  PUBLIC :: stack_backInt
  PUBLIC :: stack_backSngl
  PUBLIC :: stack_backDble
  PUBLIC :: stack_popbackInt
  PUBLIC :: stack_popbackSngl
  PUBLIC :: stack_popbackDble

  INTERFACE stack_pushback
    MODULE PROCEDURE stack_pushbackInt
    MODULE PROCEDURE stack_pushbackSngl
    MODULE PROCEDURE stack_pushbackDble
  END INTERFACE

!<types>

!<typeblock>

  ! Type block for holding a dynamic stack
  TYPE t_stack
    ! Size of stack
    INTEGER :: istackSize = 0

    ! Position of last stack item
    INTEGER :: istackPosition = 0

    ! Data type stored in stack
    INTEGER :: idataType

    ! Handle to stack data
    INTEGER :: h_StackData = ST_NOHANDLE
  END TYPE t_stack

!</typeblock>

!</types>
  
CONTAINS

!<subroutine>

  SUBROUTINE stack_create(rstack,isize,cdataType)

!<description>
    ! Creates a stack with initial memory
!</description>

!<input>
    ! Initial stack size
    INTEGER, INTENT(IN) :: isize

    ! Stack data type
    INTEGER, INTENT(IN) :: cdataType
!</input>

!<inputoutput>
    ! Stack
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>
!</subroutine>
    
    SELECT CASE (cdataType)
    CASE (ST_INT,ST_SINGLE,ST_DOUBLE)
      CALL storage_new("stack_create","h_StackData",isize,ST_DOUBLE&
          &,rstack%h_StackData,ST_NEWBLOCK_NOINIT)

    CASE DEFAULT
      PRINT *, "stack_create: Invalid data type!"
      STOP
    END SELECT

    rstack%istackSize=isize
    rstack%idataType=cdataType
  END SUBROUTINE stack_create
  
  !************************************************************************

!<subroutine>

  SUBROUTINE stack_release(rstack)

!<description>
    ! Release a stack
!</description>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>
!</subroutine>
    
    IF (rstack%h_StackData /= ST_NOHANDLE) CALL storage_free(rstack&
        &%h_StackData)
    
    rstack%istackSize=0
    rstack%istackPosition=0
  END SUBROUTINE stack_release
  
  !************************************************************************

!<subroutine>

  SUBROUTINE stack_clear(rstack)

!<description>
    ! Clear stack, i.e., reset stack pointer to zero
!</description>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>
!</subroutine>
    
    rstack%istackPosition=0       
  END SUBROUTINE stack_clear

  !************************************************************************

!<function>

  FUNCTION stack_isempty(rstack) RESULT(bisempty)

!<description>
    ! Returns TRUE if stack is empty
!</description>

!<input>
    TYPE(t_stack), INTENT(IN) :: rstack
!</input>

!<result>
    LOGICAL :: bisempty
!</result>
!</function>

    bisempty=(rstack%istackPosition==0)
  END FUNCTION stack_isempty

  !************************************************************************

!<function>

  FUNCTION stack_size(rstack) RESULT(isize)

!<description>
    ! Returns the stack size
!</description>

!<input>
    TYPE(t_stack), INTENT(IN) :: rstack
!</input>

!<result>
    INTEGER :: isize
!</result>
!</function>

    isize=rstack%istackPosition
  END FUNCTION stack_size

  !************************************************************************

!<subroutine>
  
  SUBROUTINE stack_pushbackInt(rstack,idata)

!<description>
    ! Add an integer value to the top of the stack
!</description>
    
!<input>
    INTEGER, INTENT(IN) :: idata
!</input>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: StackData

    IF (rstack%h_StackData == ST_NOHANDLE) THEN
      PRINT *, "stack_pushbackInt: Invalid data type"
      STOP
    END IF

    ! Double storage for stack if required
    IF (rstack%istackSize == rstack%istackPosition) THEN
      CALL storage_realloc("stack_pushbackInt",2*rstack%istackSize&
          &,rstack%h_StackData,ST_NEWBLOCK_NOINIT,.TRUE.)
      rstack%istackSize=2*rstack%istackSize
    END IF

    ! Increase stack pointer
    rstack%istackPosition=rstack%istackPosition+1
    
    ! Put data
    CALL storage_getbase_double(rstack%h_StackData,StackData)
    StackData(rstack%istackPosition)=TRANSFER(idata,StackData(rstack&
        &%istackPosition))
  END SUBROUTINE stack_pushbackInt

  !************************************************************************

!<subroutine>
  
  SUBROUTINE stack_pushbackSngl(rstack,sdata)

!<description>
    ! Add a single value to the top of the stack
!</description>
    
!<input>
    REAL(SP), INTENT(IN) :: sdata
!</input>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: StackData

    IF (rstack%h_StackData == ST_NOHANDLE) THEN
      PRINT *, "stack_pushbackInt: Invalid data type"
      STOP
    END IF

    ! Double storage for stack if required
    IF (rstack%istackSize == rstack%istackPosition) THEN
      CALL storage_realloc("stack_pushbackInt",2*rstack%istackSize&
          &,rstack%h_StackData,ST_NEWBLOCK_NOINIT,.TRUE.)
      rstack%istackSize=2*rstack%istackSize
    END IF

    ! Increase stack pointer
    rstack%istackPosition=rstack%istackPosition+1
    
    ! Put data
    CALL storage_getbase_double(rstack%h_StackData,StackData)
    StackData(rstack%istackPosition)=TRANSFER(sdata,StackData(rstack&
        &%istackPosition))
  END SUBROUTINE stack_pushbackSngl

  !************************************************************************

!<subroutine>
  
  SUBROUTINE stack_pushbackDble(rstack,ddata)

!<description>
    ! Add a double value to the top of the stack
!</description>
    
!<input>
    REAL(DP), INTENT(IN) :: ddata
!</input>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: StackData

    IF (rstack%h_StackData == ST_NOHANDLE) THEN
      PRINT *, "stack_pushbackInt: Invalid data type"
      STOP
    END IF

    ! Double storage for stack if required
    IF (rstack%istackSize == rstack%istackPosition) THEN
      CALL storage_realloc("stack_pushbackInt",2*rstack%istackSize&
          &,rstack%h_StackData,ST_NEWBLOCK_NOINIT,.TRUE.)
      rstack%istackSize=2*rstack%istackSize
    END IF

    ! Increase stack pointer
    rstack%istackPosition=rstack%istackPosition+1
    
    ! Put data
    CALL storage_getbase_double(rstack%h_StackData,StackData)
    StackData(rstack%istackPosition)=TRANSFER(ddata,StackData(rstack&
        &%istackPosition))
  END SUBROUTINE stack_pushbackDble

  !************************************************************************

!<function>

  FUNCTION stack_backInt(rstack) RESULT(idata)

!<description>
    ! Return integer value from top of the stack
!</description>

!<input>
    TYPE(t_stack), INTENT(IN) :: rstack
!</input>

!<result>
    INTEGER :: idata
!</result>
!</function>

    REAL(DP), DIMENSION(:), POINTER :: StackData
    
    IF (.NOT.stack_isempty(rstack)) THEN
      CALL storage_getbase_double(rstack%h_StackData,StackData)
      idata=TRANSFER(StackData(rstack%istackPosition),idata)
    ELSE
      PRINT *, "stack_backInt: Stack empty!"
      STOP
    END IF
    
  END FUNCTION stack_backInt

  !************************************************************************

!<function>

  FUNCTION stack_backSngl(rstack) RESULT(sdata)

!<description>
    ! Return single value from top of the stack
!</description>

!<input>
    TYPE(t_stack), INTENT(IN) :: rstack
!</input>

!<result>
    REAL(SP) :: sdata
!</result>
!</function>

    REAL(DP), DIMENSION(:), POINTER :: StackData
    
    IF (.NOT.stack_isempty(rstack)) THEN
      CALL storage_getbase_double(rstack%h_StackData,StackData)
      sdata=TRANSFER(StackData(rstack%istackPosition),sdata)
    ELSE
      PRINT *, "stack_backSngl: Stack empty!"
      STOP
    END IF
    
  END FUNCTION stack_backSngl

  !************************************************************************

!<function>

  FUNCTION stack_backDble(rstack) RESULT(ddata)

!<description>
    ! Return a single value from top of the stack
!</description>

!<input>
    TYPE(t_stack), INTENT(IN) :: rstack
!</input>

!<result>
    REAL(SP) :: ddata
!</result>
!</function>

    REAL(DP), DIMENSION(:), POINTER :: StackData
    
    IF (.NOT.stack_isempty(rstack)) THEN
      CALL storage_getbase_double(rstack%h_StackData,StackData)
      ddata=TRANSFER(StackData(rstack%istackPosition),ddata)
    ELSE
      PRINT *, "stack_backDble: Stack empty!"
      STOP
    END IF
    
  END FUNCTION stack_backDble

  !************************************************************************

!<function>

  FUNCTION stack_popbackInt(rstack) RESULT(idata)

!<description>
    ! Remove an integer value from top of the stack
!</description>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>

!<result>
    INTEGER :: idata
!</result>
!</function>

    REAL(DP), DIMENSION(:), POINTER :: StackData
    
    IF (.NOT.stack_isempty(rstack)) THEN
      CALL storage_getbase_double(rstack%h_StackData,StackData)
      idata=TRANSFER(StackData(rstack%istackPosition),idata)
      rstack%istackPosition=rstack%istackPosition-1
    ELSE
      PRINT *, "stack_popbackInt: Stack empty!"
      STOP
    END IF
    
  END FUNCTION stack_popbackInt

  !************************************************************************

!<function>

  FUNCTION stack_popbackSngl(rstack) RESULT(sdata)

!<description>
    ! Remove a single value from top of the stack
!</description>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>

!<result>
    REAL(SP) :: sdata
!</result>
!</function>

    REAL(DP), DIMENSION(:), POINTER :: StackData
    
    IF (.NOT.stack_isempty(rstack)) THEN
      CALL storage_getbase_double(rstack%h_StackData,StackData)
      sdata=TRANSFER(StackData(rstack%istackPosition),sdata)
      rstack%istackPosition=rstack%istackPosition-1
    ELSE
      PRINT *, "stack_popbackSngl: Stack empty!"
      STOP
    END IF
    
  END FUNCTION stack_popbackSngl

  !************************************************************************

!<function>

  FUNCTION stack_popbackDble(rstack) RESULT(ddata)

!<description>
    ! Remove a single value from top of the stack
!</description>

!<inputoutput>
    TYPE(t_stack), INTENT(INOUT) :: rstack
!</inputoutput>

!<result>
    REAL(SP) :: ddata
!</result>
!</function>

    REAL(DP), DIMENSION(:), POINTER :: StackData
    
    IF (.NOT.stack_isempty(rstack)) THEN
      CALL storage_getbase_double(rstack%h_StackData,StackData)
      ddata=TRANSFER(StackData(rstack%istackPosition),ddata)
      rstack%istackPosition=rstack%istackPosition-1
    ELSE
      PRINT *, "stack_popbackDble: Stack empty!"
      STOP
    END IF
    
  END FUNCTION stack_popbackDble
END MODULE stack
