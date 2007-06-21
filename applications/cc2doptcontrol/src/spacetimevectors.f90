!##############################################################################
!# ****************************************************************************
!# <name> spacetimevectors </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises 'space-time-vectors'. A 'space-time-vector' is
!# a list of block vectors. Every block vector represents a timestep in
!# a nonstationary simulation where the timestep length and the number of
!# timesteps is fixed.
!#
!# The module provides automatic functionality to write out vectors to
!# external memory storage in case there is too much memory allocated.
!#
!# The module provides the following subroutines:
!#
!# 1.) sptivec_initVector
!#     -> Initialise a space-time vector for a given time discretisation
!#        configuration
!#
!# 2.) sptivec_releaseVector
!#     -> Releases a space time vector.
!#
!# 3.) sptivec_setTimestepData
!#     -> Stores a timestep vector in a global space-time vector
!# 
!# 4.) sptivec_getTimestepData
!#     -> Restores a timestep vector from a global space-time vector
!#
!# 5.) sptivec_convertSupervectorToVector
!#     -> Converts a global space-time vector to a usual block vector.
!#
!# 6.) sptivec_vectorLinearComb
!#     -> Linear combination of two vectors
!#
!# 7.) sptivec_copyVector
!#     -> Copy a vector to another
!#
!# </purpose>
!##############################################################################

MODULE spacetimevectors

  USE fsystem
  USE genoutput
  USE storage
  USE linearsystemscalar
  USE linearsystemblock

  IMPLICIT NONE

!<types>

!<typeblock>

  ! This structure saves a 'space-time-vector'. A space-time-vector is basically
  ! a list of block vectors for every timestep of a nonstationary simulation
  ! which simulates simultaneously in space and time. Some parts of this
  ! vector may be written to disc if there's not enough memory.
  TYPE t_spacetimeVector
  
    ! The name of a directory that is used for temporarily storing data to disc.
    ! Standard is the current directory.
    CHARACTER(SYS_STRLEN) :: sdirectory = './'
    
    ! Number of equations in each subvector of the 'global time-step vector'.
    INTEGER(PREC_VECIDX) :: NEQ = 0
    
    ! Number of subvectors/timesteps saved in p_IdataHandleList.
    INTEGER :: ntimesteps
    
    ! A list of handles (dimension 0:ntimesteps) to double-precision arrays 
    ! which save the data of the ntimesteps+1 data subvectors. 
    ! A handle number > 0 indicates that the vector data
    ! is appearent in memory. A handle number < 0 indicates that the array
    ! is written to disc and has to be read in before it can be used.
    ! A value ST_NOHANDLE indicates that there is currently no data stored
    ! at that timestep.
    INTEGER, DIMENSION(:), POINTER :: p_IdataHandleList => NULL()
    
  END TYPE

!</typeblock>

!</types>

CONTAINS

!<subroutine>

  SUBROUTINE sptivec_initVector (rspaceTimeVector,NEQ,ntimesteps,sdirectory)

!<description>
  ! Initialises a space time vector. rvectorTemplate is a template vector that
  ! defines the basic shape of the data vectors in all time steps that are
  ! to be maintained. ntimesteps defines the number of timesteps to maintain.
!</desctiprion>

!<input>
  ! Number of equations in the vectors
  INTEGER(PREC_VECIDX), INTENT(IN) :: NEQ
  
  ! Number of timesteps to maintain.
  ! The number of subvectors that is reserved is therefore ntimesteps+1!
  INTEGER, INTENT(IN) :: ntimesteps
  
  ! OPTIONAL: Directory name where to save temporary data. If not present,
  ! the current directory './' is the default.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sdirectory
!</input>

!<output>
  ! Space-time vector structure to be initialised.
  TYPE(t_spacetimeVector), INTENT(OUT) :: rspaceTimeVector
!</output>

!</subroutine>

    INTEGER :: i

    ! Initialise the data.
    IF (PRESENT(sdirectory)) rspaceTimeVector%sdirectory = sdirectory
    rspaceTimeVector%ntimesteps = ntimesteps
    ALLOCATE(rspaceTimeVector%p_IdataHandleList(0:ntimesteps))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    
    rspaceTimeVector%NEQ = NEQ
    
    ! Allocate memory for every subvector
    DO i=0,ntimesteps
      CALL storage_new ('sptivec_initVector', 'stvec_'//TRIM(sys_si(i,10)), &
        NEQ, ST_DOUBLE, rspaceTimeVector%p_IdataHandleList(i), ST_NEWBLOCK_ZERO)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_releaseVector (rspaceTimeVector)

!<description>
  ! Releases a space-time vector. All allocated memory is released. Temporary data
  ! on disc is deleted.
!</desctiprion>

!<inputoutput>
  ! Space-time vector structure to be initialised.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rspaceTimeVector
!</inputoutput>

!</subroutine>

    ! local variables -- initialised by Fortran default initialisation!
    TYPE(t_spacetimeVector) :: rspaceTimeVectorTempl
    INTEGER :: i

    ! Deallocate data
    DO i=0,rspaceTimeVector%ntimesteps
      IF (rspaceTimeVector%p_IdataHandleList(i) .NE. ST_NOHANDLE) THEN
      
        IF (rspaceTimeVector%p_IdataHandleList(i) .GT. 0) THEN
          CALL storage_free (rspaceTimeVector%p_IdataHandleList(i))
        ELSE IF (rspaceTimeVector%p_IdataHandleList(i) .LT. 0) THEN
          PRINT *,'external data not implemented!'
          STOP
        END IF
        
      END IF
    END DO

    DEALLOCATE(rspaceTimeVector%p_IdataHandleList)

    ! Initialise with default values.
    rspaceTimeVector = rspaceTimeVectorTempl

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_setTimestepData (rspaceTimeVector, isubvector, rvector)

!<description>
  ! Stores the data of rvector at timestep itimestep in the rspaceTimeVector.
  ! structure.
!</desctiprion>

!<input>
  ! Number of the subvector that corresponds to rvector. >= 0, <= ntimesteps from
  ! the initialisation of rspaceTimeVector.
  INTEGER, INTENT(IN) :: isubvector
  
  ! Vector with data that should be associated to timestep itimestep.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
!</input>

!<inputoutput>
  ! Space-time vector structure where to save the data.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rspaceTimeVector
!</inputoutput>

!</subroutine>

    ! Make sure we can store the timestep data.
    IF ((isubvector .LT. 0) .OR. (isubvector .GT. rspaceTimeVector%ntimesteps)) THEN
      PRINT *,'Invalid timestep number!'
      STOP
    END IF
    
    IF (rvector%NEQ .NE. rspaceTimeVector%NEQ) THEN
      PRINT *,'Vector size invalid!'
      STOP
    END IF

    IF ((rspaceTimeVector%p_IdataHandleList(isubvector) .NE. ST_NOHANDLE) .AND. &
        (rspaceTimeVector%p_IdataHandleList(isubvector) .LT. 0)) THEN
      PRINT *,'external data not implemented!'
      STOP
    END IF

    ! Save the vector data. If necessary, new memory is allocated -- as the
    ! default value of the handle in the handle list is ST_NOHANDLE.
    CALL storage_copy (rvector%h_Ddata,rspaceTimeVector%p_IdataHandleList(isubvector))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_getTimestepData (rspaceTimeVector, isubvector, rvector)

!<description>
  ! Restores the data of timestep itimestep into the vector rvector.
!</desctiprion>

!<input>
  ! Number of the subvector that corresponds to rvector. >= 0, <= ntimesteps from
  ! the initialisation of rspaceTimeVector.
  INTEGER, INTENT(IN) :: isubvector
  
  ! Space-time vector structure where to save the data.
  TYPE(t_spacetimeVector), INTENT(IN) :: rspaceTimeVector
!</input>

!<inputoutput>
  ! Vector with data that should receive the data.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! Make sure we can store the timestep data.
    IF ((isubvector .LT. 0) .OR. (isubvector .GT. rspaceTimeVector%ntimesteps)) THEN
      PRINT *,'Invalid timestep number!'
      STOP
    END IF
    
    IF (rvector%NEQ .NE. rspaceTimeVector%NEQ) THEN
      PRINT *,'Vector size invalid!'
      STOP
    END IF
    
    IF ((rspaceTimeVector%p_IdataHandleList(isubvector) .NE. ST_NOHANDLE) .AND. &
        (rspaceTimeVector%p_IdataHandleList(isubvector) .LT. 0)) THEN
      PRINT *,'external data not implemented!'
      STOP
    END IF

    ! Save the vector data. If necessary, new memory is allocated -- as the
    ! default value of the handle in the handle list is ST_NOHANDLE.
    CALL storage_copy (rspaceTimeVector%p_IdataHandleList(isubvector),rvector%h_Ddata)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_vectorLinearComb (rx,ry,cx,cy)

!<description>
  ! Performs a linear combination of space-time vectors: ry = cx * rx  +  cy * ry
!</desctiprion>

!<input>
  ! First source vector
  TYPE(t_spacetimeVector), INTENT(IN)   :: rx
  
  ! Scaling factor for Dx
  REAL(DP), INTENT(IN)               :: cx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)               :: cy
!</input>

!<inputoutput>
  ! Second source vector; also receives the result
  TYPE(t_spacetimeVector), INTENT(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

    INTEGER :: i
    INTEGER(PREC_VECIDX), DIMENSION(1) :: Isize
    TYPE(t_vectorBlock) :: rxBlock,ryBlock
    
    IF (rx%NEQ .NE. ry%NEQ) THEN
      PRINT *,'Space-time vectors have different size!'
      STOP
    END IF

    IF (rx%ntimesteps .NE. ry%ntimesteps) THEN
      PRINT *,'Space-time vectors have different number of timesteps!'
      STOP
    END IF
    
    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    CALL lsysbl_createVecBlockDirect (rxBlock,Isize,.FALSE.)
    CALL lsysbl_createVecBlockDirect (ryBlock,Isize,.FALSE.)

    ! Loop through the substeps, load the data in, perform the linear combination
    ! and write out again.
    DO i=0,rx%ntimesteps
      CALL sptivec_getTimestepData (rx, i, rxBlock)
      CALL sptivec_getTimestepData (ry, i, ryBlock)
      
      CALL lsysbl_vectorLinearComb (rxBlock,ryBlock,cx,cy)

      CALL sptivec_setTimestepData (ry, i, ryBlock)
    END DO

    ! Release temp memory    
    CALL lsysbl_releaseVector (ryBlock)
    CALL lsysbl_releaseVector (rxBlock)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_copyVector (rx,ry)

!<description>
  ! Copys a vector: ry := rx.
!</desctiprion>

!<input>
  ! Source vector
  TYPE(t_spacetimeVector), INTENT(IN)   :: rx
!</input>

!<inputoutput>
  ! Destination vector
  TYPE(t_spacetimeVector), INTENT(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

    INTEGER :: i
    INTEGER(PREC_VECIDX), DIMENSION(1) :: Isize
    TYPE(t_vectorBlock) :: rxBlock,ryBlock
    
    IF (rx%NEQ .NE. ry%NEQ) THEN
      PRINT *,'Space-time vectors have different size!'
      STOP
    END IF

    IF (rx%ntimesteps .NE. ry%ntimesteps) THEN
      PRINT *,'Space-time vectors have different number of timesteps!'
      STOP
    END IF
    
    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    CALL lsysbl_createVecBlockDirect (rxBlock,Isize,.FALSE.)
    CALL lsysbl_createVecBlockDirect (ryBlock,Isize,.FALSE.)

    ! Loop through the substeps, load the data in, perform the linear combination
    ! and write out again.
    DO i=1,rx%ntimesteps
      CALL sptivec_getTimestepData (rx, i, rxBlock)
      CALL sptivec_getTimestepData (ry, i, ryBlock)
      
      CALL lsysbl_copyVector (rxBlock,ryBlock)

      CALL sptivec_setTimestepData (ry, i, ryBlock)
    END DO

    ! Release temp memory    
    CALL lsysbl_releaseVector (ryBlock)
    CALL lsysbl_releaseVector (rxBlock)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_convertSupervecToVector (rxsuper, rx)

!<description>
  ! Converts a space-time coupled vector into a usual block vector.
  ! The blocks in rx correspond to the timesteps in rxsuper.
  !
  ! WARNING: This may be memory intensive!
!</description>

!<input>
  ! Space time vector to be converted.
  TYPE(t_spacetimeVector), INTENT(IN)   :: rxsuper
!</input>

!<inputoutput>
  ! Destination block vector that should receive the result.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rx
!</inputoutput>

!</subroutine>

    INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: Isize
    TYPE(t_vectorBlock) :: rvectorTmp
    INTEGER :: i
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata1,p_Ddata2
    
    ! Create a complete new vector
    IF (rx%NEQ .NE. 0) CALL lsysbl_releaseVector (rx)
  
    ! Create a vector in the correct size
    ALLOCATE(Isize(rxsuper%ntimesteps+1))
    Isize(:) = rxsuper%NEQ
    CALL lsysbl_createVecBlockDirect (rx,Isize,.FALSE.)

    ! Create a 1-block temp vector for the data    
    CALL lsysbl_createVecBlockDirect (rvectorTmp,Isize(1:1),.FALSE.)
    CALL lsysbl_getbase_double (rvectorTmp,p_Ddata1)
    
    ! Load the subvectors and write them to the global vector.
    DO i=0,rxsuper%ntimesteps
      CALL sptivec_getTimestepData (rxsuper, i, rvectorTmp)
      
      CALL lsyssc_getbase_double (rx%RvectorBlock(i+1),p_Ddata2)
      CALL lalg_copyVectorDble (p_Ddata1,p_Ddata2)
    END DO
    
    ! Release the temp vector
    CALL lsysbl_releaseVector (rvectorTmp)
    
    DEALLOCATE(Isize)
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_convertVectorToSupervec (rx,rxsuper)

!<description>
  ! Converts a usual block vector into a space-time coupled vector.
  ! The blocks in rx correspond to the timesteps in rxsuper.
!</description>

!<input>
  ! Block vector to be converted.
  TYPE(t_vectorBlock), INTENT(IN)   :: rx
!</input>

!<inputoutput>
  ! Destination space-time vector that should receive the result.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rxsuper
!</inputoutput>

!</subroutine>

    INTEGER(PREC_VECIDX), DIMENSION(1) :: Isize
    TYPE(t_vectorBlock) :: rvectorTmp
    INTEGER :: i
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata1,p_Ddata2
    
    ! Create a 1-block temp vector for the data    
    Isize(1) = rxsuper%NEQ
    CALL lsysbl_createVecBlockDirect (rvectorTmp,Isize(1:1),.FALSE.)
    CALL lsysbl_getbase_double (rvectorTmp,p_Ddata2)
    
    ! Load the subvectors and write them to the global vector.
    DO i=0,rxsuper%ntimesteps
      CALL lsyssc_getbase_double (rx%RvectorBlock(i+1),p_Ddata1)
      CALL lalg_copyVectorDble (p_Ddata1,p_Ddata2)

      CALL sptivec_setTimestepData (rxsuper, i, rvectorTmp)
    END DO
    
    ! Release the temp vector
    CALL lsysbl_releaseVector (rvectorTmp)
    
  END SUBROUTINE

END MODULE
