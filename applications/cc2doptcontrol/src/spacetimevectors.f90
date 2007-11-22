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
!# 5.) sptivec_getTimestepDataByTime
!#     -> Restores a timestep vector from a global space-time vector based
!#        on a time stamp
!#
!# 6.) sptivec_convertSupervectorToVector
!#     -> Converts a global space-time vector to a usual block vector.
!#
!# 7.) sptivec_vectorLinearComb
!#     -> Linear combination of two vectors
!#
!# 8.) sptivec_copyVector
!#     -> Copy a vector to another
!#
!# 9.) sptivec_vectorNorm
!#     -> Calculate the norm of a vector
!#
!# 10.) sptivec_clearVector
!#      -> Clears a vector
!#
!# 11.) sptivec_saveToCollection
!#      -> Saves a space-time vector to a collection
!#
!# 12.) sptivec_restoreFromCollection
!#      -> Restores a space-time vector from a collection
!#
!# 13.) sptivec_removeFromCollection
!#      -> Removes a space-time vector from a collection
!#
!# 14.) sptivec_setConstant
!#      -> Initialises the whole space-time vector with a constant value.
!#
!# 15.) sptivec_loadFromFileSequence
!#      -> Reads in a space-time vector from a sequence of files on disc
!#
!# 16.) sptivec_saveToFileSequence
!#      -> Writes a space-time vector to a sequence of files on disc
!#
!# </purpose>
!##############################################################################

MODULE spacetimevectors

  USE fsystem
  USE genoutput
  USE storage
  USE linearsystemscalar
  USE linearsystemblock
  USE collection
  USE vectorio

  IMPLICIT NONE

!<types>

!<typeblock>

  ! This structure saves a 'space-time-vector'. A space-time-vector is basically
  ! a list of block vectors for every timestep of a nonstationary simulation
  ! which simulates simultaneously in space and time. Some parts of this
  ! vector may be written to disc if there's not enough memory.
  !
  ! The vector can be accessed in two ways: On the one hand, one can read a
  ! specific time step by its number. On the other hand, one can access the
  ! content by a time stamp.
  TYPE t_spacetimeVector
  
    ! Whether this vector shares its data with another vector.
    LOGICAL :: bisCopy = .FALSE.
  
    ! The name of a directory that is used for temporarily storing data to disc.
    ! Standard is the current directory.
    CHARACTER(SYS_STRLEN) :: sdirectory = './'
    
    ! This flag defines a scaling factor for each substep. The scaling is applied
    ! when a subvector is read and is reset to 1.0 if a subvector is saved.
    ! The whole vector is zero if it's new and if it's cleared by clearVector.
    REAL(DP), DIMENSION(:), POINTER :: p_Dscale => NULL()
    
    ! Number of equations in each subvector of the 'global time-step vector'.
    INTEGER(PREC_VECIDX) :: NEQ = 0
    
    ! Number of subvectors/timesteps saved in p_IdataHandleList.
    INTEGER :: ntimesteps = 0
    
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

  ! ***************************************************************************

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
    ALLOCATE(rspaceTimeVector%p_Dscale(0:ntimesteps))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    rspaceTimeVector%p_Dscale(:) = 0.0_DP
    
    rspaceTimeVector%NEQ = NEQ
    
    ! Allocate memory for every subvector
    DO i=0,ntimesteps
      CALL storage_new ('sptivec_initVector', 'stvec_'//TRIM(sys_siL(i,10)), &
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

    IF (.NOT. ASSOCIATED(rspaceTimeVector%p_IdataHandleList)) THEN
      CALL output_line('Warning: Releasing unused vector!',&
          ssubroutine='sptivec_releaseVector')
      RETURN
    END IF

    ! Deallocate data -- if the data is not shared with another vector
    IF (.NOT. rspaceTimeVector%bisCopy) THEN
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
    END IF

    DEALLOCATE(rspaceTimeVector%p_IdataHandleList)
    DEALLOCATE(rspaceTimeVector%p_Dscale)

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
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest

    ! Make sure we can store the timestep data.
    IF ((isubvector .LT. 0) .OR. (isubvector .GT. rspaceTimeVector%ntimesteps)) THEN
      CALL output_line('Invalid timestep number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_setTimestepData')
      CALL sys_halt()
    END IF
    
    IF (rvector%NEQ .NE. rspaceTimeVector%NEQ) THEN
      CALL output_line('Vector size invalid!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_setTimestepData')
      CALL sys_halt()
    END IF

    IF ((rspaceTimeVector%p_IdataHandleList(isubvector) .NE. ST_NOHANDLE) .AND. &
        (rspaceTimeVector%p_IdataHandleList(isubvector) .LT. 0)) THEN
      CALL output_line('external data not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_setTimestepData')
      CALL sys_halt()
    END IF

    ! Save the vector data. If necessary, new memory is allocated -- as the
    ! default value of the handle in the handle list is ST_NOHANDLE.
    !
    ! Don't use storage_copy, as this might give errors in case the array
    ! behind the handle is longer than the vector!
    CALL lsysbl_getbase_double (rvector,p_Dsource)
    CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(isubvector),p_Ddest)
    CALL lalg_copyVectorDble (p_Dsource,p_Ddest)

    ! After a setTimestepData, the scale factor is 1.0.
    rspaceTimeVector%p_Dscale(isubvector) = 1.0_DP

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

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest

    ! Make sure we can store the timestep data.
    IF ((isubvector .LT. 0) .OR. (isubvector .GT. rspaceTimeVector%ntimesteps)) THEN
      CALL output_line('Invalid timestep number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      CALL sys_halt()
    END IF
    
    IF (rvector%NEQ .NE. rspaceTimeVector%NEQ) THEN
      CALL output_line('Vector size invalid!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      CALL sys_halt()
    END IF
    
    IF (rspaceTimeVector%p_Dscale(isubvector) .EQ. 0.0_DP) THEN
     ! The vector is a zero vector
      CALL lsysbl_clearVector (rvector)
      RETURN
    END IF

    IF ((rspaceTimeVector%p_IdataHandleList(isubvector) .NE. ST_NOHANDLE) .AND. &
        (rspaceTimeVector%p_IdataHandleList(isubvector) .LT. 0)) THEN
      CALL output_line('external data not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      CALL sys_halt()
    END IF
    
    ! Get the vector data. 
    !
    ! Don't use storage_copy, as this might give errors in case the array
    ! behind the handle is longer than the vector!
    CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(isubvector),&
        p_Dsource)
    CALL lsysbl_getbase_double (rvector,p_Ddest)
    CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
    
    ! Scale the vector?
    IF (rspaceTimeVector%p_Dscale(isubvector) .NE. 1.0_DP) THEN
      CALL lalg_scaleVectorDble (p_Ddest,rspaceTimeVector%p_Dscale(isubvector))
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_getTimestepDataByTime (rspaceTimeVector, dtimestamp, rvector)

!<description>
  ! Restores the data of a timestep into the vector rvector. dtimestamp is a time
  ! stamp in the range 0.0 .. 1.0, where 0.0 corresponds to the 0th subvector
  ! and 1.0 to the last subvector in the space time vector. If dtimestamp
  ! specifies a time stamp between two stored vectors, quadratic interpolation
  ! is used to calculate rvector.
!</desctiprion>

!<input>
  ! Time stamp of the vector whose data should be retrieved.
  REAL(DP), INTENT(IN) :: dtimestamp
  
  ! Space-time vector structure where to save the data.
  TYPE(t_spacetimeVector), INTENT(IN) :: rspaceTimeVector
!</input>

!<inputoutput>
  ! Vector with data that should receive the data.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Dsource1,p_Dsource2,p_Dsource3,p_Ddest
    INTEGER :: itimestep1,itimestep2,itimestep3
    REAL(DP) :: dreltime,dabstime
    INTEGER :: i
    REAL(DP) :: dscal1,dscal2,dscal3

    ! Make sure we can store the timestep data.
    IF ((dtimestamp .LT. 0.0_DP) .OR. (dtimestamp .GT. 1.0_DP)) THEN
      CALL output_line('Invalid time stamp!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepDataByTime')
      CALL sys_halt()
    END IF
    
    IF (rvector%NEQ .NE. rspaceTimeVector%NEQ) THEN
      CALL output_line('Vector size invalid!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      CALL sys_halt()
    END IF
    
    ! Get the time step which is closest to the time stamp.
    ! Rescale dtimestamp to the interval [0..ntimesteps].
    dabstime = dtimestamp*REAL(rspaceTimeVector%ntimesteps,DP)
    itimestep2 = INT(dabstime + 0.5_DP)
    
    IF (dabstime .EQ. REAL(itimestep2,DP)) THEN
      ! Nice coincidence, we have exactly timestep itimestep2. Ok, then we 
      ! can call the routine to get that timestep; this saves us some
      ! time as the interpolation can be omitted.
      CALL sptivec_getTimestepData (rspaceTimeVector, itimestep2, rvector)
      RETURN
    END IF
    
    IF (rspaceTimeVector%ntimesteps .EQ. 1) THEN
      ! Special case: only one timestep!
      itimestep1 = 0
      itimestep2 = 0
      itimestep3 = 1
    ELSE
      ! Is this the first or the last timestep?
      IF (itimestep2 .EQ. 0) THEN
        ! First timestep. Interpolate between timesteps 0,1 and 2, evaluate 
        ! near timestep 0.
        itimestep1 = 0
        itimestep2 = 1
        itimestep3 = 2
      ELSE IF (itimestep2 .EQ. rspaceTimeVector%ntimesteps) THEN
        ! Last timestep. Interpolate between timesteps n-2,n-1 and n, evaluate 
        ! near timestep n.
        itimestep1 = rspaceTimeVector%ntimesteps-2
        itimestep2 = rspaceTimeVector%ntimesteps-1
        itimestep3 = rspaceTimeVector%ntimesteps
      ELSE
        ! Somewhere in the inner. Get the number of the previous and next timestep
        itimestep1 = itimestep2-1
        itimestep3 = itimestep2+1
      END IF
    END IF

    ! Calculate the 'relative' time in the interval [-1,1], where -1 corresponds
    ! to timestep itimestep1, 0 to itimestep2 and +1 to itimestep3.
    ! This will be used to evaluate the quadratic polynomial.
    dreltime = dabstime-REAL(itimestep2,DP)
    
    IF (((rspaceTimeVector%p_IdataHandleList(itimestep1) .NE. ST_NOHANDLE) .AND. &
         (rspaceTimeVector%p_IdataHandleList(itimestep1) .LT. 0)) .OR. &
        ((rspaceTimeVector%p_IdataHandleList(itimestep2) .NE. ST_NOHANDLE) .AND. &
         (rspaceTimeVector%p_IdataHandleList(itimestep2) .LT. 0)) .OR. &
        ((rspaceTimeVector%p_IdataHandleList(itimestep3) .NE. ST_NOHANDLE) .AND. &
         (rspaceTimeVector%p_IdataHandleList(itimestep3) .LT. 0))) THEN
      CALL output_line('external data not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      CALL sys_halt()
    END IF

    ! Get the vector data of the three timesteps
    CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(itimestep1),&
        p_Dsource1)
    CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(itimestep2),&
        p_Dsource2)
    CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(itimestep3),&
        p_Dsource3)

    CALL lsysbl_getbase_double (rvector,p_Ddest)

    ! Calculate the quadratic interpolation of the three arrays.
    dscal1 = rspaceTimeVector%p_Dscale(itimestep1)
    dscal2 = rspaceTimeVector%p_Dscale(itimestep2)
    dscal3 = rspaceTimeVector%p_Dscale(itimestep3)
    IF (rspaceTimeVector%ntimesteps .EQ. 1) THEN
      ! Special case: only 1 timestep. Linear interpolation. dreltime is in [0..1]!
      DO i=1,SIZE(p_Ddest)
        p_Ddest(i) = (1.0_DP-dreltime) * dscal2*p_Dsource2(i) + &
                     dreltime*dscal3*p_Dsource3(i)
      END DO
    ELSE
      ! Quadratic interpolation
      DO i=1,SIZE(p_Ddest)
        CALL mprim_quadraticInterpolation (dreltime,&
            dscal1*p_Dsource1(i),dscal2*p_Dsource2(i),dscal3*p_Dsource3(i),p_Ddest(i))
      END DO
    END IF
    
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

    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dy
    
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

    ! DEBUG!!!
    CALL lsysbl_getbase_double (rxBlock,p_Dx)
    CALL lsysbl_getbase_double (ryBlock,p_Dy)

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
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dy
    
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
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rxBlock,p_Dx)
    CALL lsysbl_getbase_double (ryBlock,p_Dy)

    ! Loop through the substeps, load the data in, perform the linear combination
    ! and write out again.
    DO i=0,rx%ntimesteps
      
      IF (rx%p_Dscale(i) .EQ. 0.0_DP) THEN
        ry%p_Dscale(i) = 0.0_DP
      ELSE
        CALL sptivec_getTimestepData (rx, i, rxBlock)

        CALL lsysbl_copyVector (rxBlock,ryBlock)

        CALL sptivec_setTimestepData (ry, i, ryBlock)
      END IF

    END DO

    ! Release temp memory    
    CALL lsysbl_releaseVector (ryBlock)
    CALL lsysbl_releaseVector (rxBlock)
      
  END SUBROUTINE

  ! ***************************************************************************

!<function>
  
  REAL(DP) FUNCTION sptivec_scalarProduct (rx, ry)
  
!<description>
  ! Calculates a scalar product of two block vectors.
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>

!<input>
  ! First source vector
  TYPE(t_spacetimeVector), INTENT(IN)   :: rx

  ! Second source vector
  TYPE(t_spacetimeVector), INTENT(IN)   :: ry
!</input>

!<result>
  ! The scalar product (rx,ry) of the two block vectors.
!</result>

!</function>

    INTEGER :: i
    INTEGER(PREC_VECIDX), DIMENSION(1) :: Isize
    REAL(DP) :: dres
    TYPE(t_vectorBlock) :: rxBlock,ryBlock
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dy
    
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
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rxBlock,p_Dx)
    CALL lsysbl_getbase_double (ryBlock,p_Dy)

    ! Loop through the substeps, load the data in, perform the scalar product.
    dres = 0.0_DP
    DO i=0,rx%ntimesteps
      
      IF ((rx%p_Dscale(i) .NE. 0.0_DP) .AND. (ry%p_Dscale(i) .NE. 0.0_DP)) THEN
        CALL sptivec_getTimestepData (rx, i, rxBlock)
        CALL sptivec_getTimestepData (ry, i, ryBlock)
        
        dres = dres + lsysbl_scalarProduct (rxBlock,ryBlock)
      END IF

    END DO

    ! Release temp memory    
    CALL lsysbl_releaseVector (ryBlock)
    CALL lsysbl_releaseVector (rxBlock)
      
    sptivec_scalarProduct = dres
      
  END FUNCTION

  ! ***************************************************************************

!<function>

  REAL(DP) FUNCTION sptivec_vectorNorm (rx,cnorm)

!<description>
  ! Calculates the norm of the vector rx.
!</desctiprion>

!<input>
  ! Source vector
  TYPE(t_spacetimeVector), INTENT(IN)   :: rx

  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  INTEGER, INTENT(IN) :: cnorm
!</input>

!<result>
  ! Norm of the vector.
!</result>
  
!</function>

    INTEGER(PREC_VECIDX), DIMENSION(1) :: Isize
    TYPE(t_vectorBlock) :: rxBlock
    REAL(DP) :: dnorm
    INTEGER :: i

    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    CALL lsysbl_createVecBlockDirect (rxBlock,Isize,.FALSE.)
    
    dnorm = 0.0_DP
    
    ! Loop through the substeps, load the data in, sum up to the norm.
    DO i=0,rx%ntimesteps
      IF (rx%p_Dscale(i) .NE. 0.0_DP) THEN
        CALL sptivec_getTimestepData (rx, i, rxBlock)
        
        SELECT CASE (cnorm)
        CASE (LINALG_NORML2)
          dnorm = dnorm + lsysbl_vectorNorm (rxBlock,cnorm)**2
        CASE DEFAULT
          dnorm = dnorm + lsysbl_vectorNorm (rxBlock,cnorm)
        END SELECT
      END IF
    END DO

    ! Release temp memory    
    CALL lsysbl_releaseVector (rxBlock)
    
    ! Calculate the actual norm.
    SELECT CASE (cnorm)
    CASE (LINALG_NORML1)
      sptivec_vectorNorm = dnorm / (rx%ntimesteps+1)
    CASE (LINALG_NORML2)
      sptivec_vectorNorm = SQRT(dnorm / (rx%ntimesteps+1))
    CASE DEFAULT
      sptivec_vectorNorm = dnorm
    END SELECT

  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_clearVector (rx)

!<description>
  ! Initialises a vector with zero.
!</desctiprion>

!<inputoutput>
  ! Vector to be cleared.
  TYPE(t_spacetimeVector), INTENT(INOUT)   :: rx
!</inputoutput>

!</subroutine>

    ! Local variables
    !!REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    !!INTEGER :: i
    !!
    ! Loop over the files
    !!DO i=0,rx%ntimesteps
    !!
    !!  ! Get the data and set to a defined value.
    !!  CALL storage_getbase_double (rx%p_IdataHandleList(i),p_Ddata)
    !!  p_Ddata(:) = 0.0_DP
    !!
    !!END DO

    ! Simply set the "empty" flag to TRUE.
    ! When restoreing data with getTimestepData, that routine will return a zero vector.
    rx%p_Dscale(:) = 0.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_scaleVector (rx,dscale)

!<description>
  ! Scales a vector by dscale.
!</desctiprion>

!<input>
  ! Scaling factor.
  REAL(DP), INTENT(IN) :: dscale
!</input>

!<inputoutput>
  ! Vector to be scaled.
  TYPE(t_spacetimeVector), INTENT(INOUT)   :: rx
!</inputoutput>

!</subroutine>

    ! Local variables
    !!REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    !!INTEGER :: i
    !!
    ! Loop over the files
    !!DO i=0,rx%ntimesteps
    !!
    !!  ! Get the data and set to a defined value.
    !!  CALL storage_getbase_double (rx%p_IdataHandleList(i),p_Ddata)
    !!  p_Ddata(:) = dscale*p_Ddata(:)
    !!
    !!END DO

    ! Scale the scaling factors of all subvectors with dscale.
    rx%p_Dscale(:) = rx%p_Dscale(:) * dscale

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

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_saveToCollection (rx,rcollection,sname,ilevel,ssection)

!<description>
  ! Saves a space-time vector to a collection.
!</description>

!<input>
  ! Space-time vector to be saved.
  TYPE(t_spacetimeVector), INTENT(IN)   :: rx
  
  ! Name that the vector should be given in the collection.
  CHARACTER(LEN=*), INTENT(IN) :: sname

  ! OPTIONAL: Name of the section in the collection where to save the vector to.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssection
  
  ! OPTIONAL: Level in the collection where to save the vector to.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel
!</input>

!<inputoutput>
  ! Collection structure. The vector is saved to this.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</inputoutput>

!</subroutine>

    ! Save the content of the structure    
    CALL collct_setvalue_string (rcollection,TRIM(sname)//'_DIR',&
        rx%sdirectory,.TRUE.,ilevel,ssection)

    CALL collct_setvalue_int (rcollection,TRIM(sname)//'_NEQ',&
        rx%NEQ,.TRUE.,ilevel,ssection)
  
    CALL collct_setvalue_int (rcollection,TRIM(sname)//'_NTST',&
        rx%ntimesteps,.TRUE.,ilevel,ssection)

    IF (rx%ntimesteps .NE. 0) THEN
      CALL collct_setvalue_intarr (rcollection,TRIM(sname)//'_NTST',&
          rx%p_IdataHandleList,.TRUE.,ilevel,ssection)
      
      CALL collct_setvalue_realarr (rcollection,TRIM(sname)//'_SCALE',&
          rx%p_Dscale,.TRUE.,ilevel,ssection)

      ! Otherwise the pointers are NULL()!
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_restoreFromCollection (rx,rcollection,sname,ilevel,ssection)

!<description>
  ! Restores a space-time vector from a collection.
  !
  ! Note that this creates a copy of a previous space-time vector. This vector
  ! must be released with releaseVector to prevent memory leaks
!</description>

!<input>
  ! Collection structure where to restore rx from.
  TYPE(t_collection), INTENT(INOUT) :: rcollection

  ! Name that the vector should be given in the collection.
  CHARACTER(LEN=*), INTENT(IN) :: sname

  ! OPTIONAL: Name of the section in the collection where to save the vector to.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssection
  
  ! OPTIONAL: Level in the collection where to save the vector to.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel
!</input>

!<inputoutput>
  ! Space-time vector that receives the data.
  TYPE(t_spacetimeVector), INTENT(OUT)   :: rx
!</inputoutput>

!</subroutine>
    
    ! The vector is a copy of another one, so the releaseVector routine
    ! will not release the content.
    rx%bisCopy = .TRUE.

    ! Get the content of the structure    
    CALL collct_getvalue_string (rcollection,TRIM(sname)//'_DIR',&
        rx%sdirectory,ilevel,ssection)

    rx%NEQ = collct_getvalue_int (rcollection,TRIM(sname)//'_NEQ',&
        ilevel,ssection)
  
    rx%ntimesteps = collct_getvalue_int (rcollection,TRIM(sname)//'_NTST',&
        ilevel,ssection)
        
    IF (rx%ntimesteps .NE. 0) THEN
      ! For the handle list, we need to allocate some memory...
      ALLOCATE(rx%p_IdataHandleList(rx%ntimesteps))

      CALL collct_getvalue_intarr (rcollection,TRIM(sname)//'_NTST',&
          rx%p_IdataHandleList,ilevel,ssection)

      CALL collct_getvalue_realarr (rcollection,TRIM(sname)//'_SCALE',&
          rx%p_Dscale,ilevel,ssection)
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_removeFromCollection (rcollection,sname,ilevel,ssection)

!<description>
  ! Removes a space-time vector from a collection.
!</description>

!<input>
  ! Name that the vector should be given in the collection.
  CHARACTER(LEN=*), INTENT(IN) :: sname

  ! OPTIONAL: Name of the section in the collection where to save the vector to.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssection
  
  ! OPTIONAL: Level in the collection where to save the vector to.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel
!</input>

!<inputoutput>
  ! Collection structure where to remove sname from.
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</inputoutput>

!</subroutine>

    ! Remove all entries
    CALL collct_deleteValue (rcollection,TRIM(sname)//'_DIR',&
        ilevel,ssection)

    CALL collct_deleteValue (rcollection,TRIM(sname)//'_SCALE',&
        ilevel,ssection)
  
    CALL collct_deleteValue (rcollection,TRIM(sname)//'_NEQ',&
        ilevel,ssection)
  
    CALL collct_deleteValue (rcollection,TRIM(sname)//'_NTST',&
        ilevel,ssection)
        
    CALL collct_deleteValue (rcollection,TRIM(sname)//'_NTST',&
        ilevel,ssection)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_loadFromFileSequence (rx,sfilename,istart,iend,idelta,&
      bformatted,brepeatLast,sdirectory)

!<description>
  ! This routine loads a space-time vector from a sequence of files on the
  ! hard disc. sfilename is a directory/file name pattern in the format of
  ! a FORMAT statement that forms the filename; this pattern must contain
  ! exactly one integer format specifier, which is replaced by the
  ! file number in this routine (e.g. ' (''vector.txt.'',I5.5) ' will
  ! load a file sequence 'vector.txt.00001','vector.txt.00002',
  ! 'vector.txt.00003', ...).
  ! istart and iend prescribe the minimum/maximum file number that is
  ! inserted into the filename: The routine loads in all all files
  ! from istart to iend and forms a space-time vector from that.
  !
  ! The source files are expected to have been written out by
  ! vecio_writeBlockVectorHR or vecio_writeVectorHR.
  !
  ! If a file in the sequence is missing, that file is skipped. The corresponding
  ! solution is set to zero.
!</description>

!<input>
  ! Filename pattern + path where to form a filename from.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename

  ! Number of the first file to be read in
  INTEGER, INTENT(IN) :: istart
  
  ! Number of the last file to be read in
  INTEGER, INTENT(IN) :: iend
  
  ! Delta parameter that specifies how to increase the filename suffix.
  ! Standard is =1.
  INTEGER, INTENT(IN) :: idelta
  
  ! Whether to read formatted or unformatted data from disc.
  LOGICAL, INTENT(IN) :: bformatted

  ! OPTIONAL: Repetition of last solution.
  ! If this value is set to TRUE and there are not enough solutions
  ! saved on disc for all timesteps, the last available solution is used
  ! to fill up the missing solutions.
  ! Can be used to continue a quasi-stationary solution if only the
  ! first couple of solutions up to a point in time are saved while
  ! the remaining (not available on disc) solutions are identical
  ! to the last one.
  !
  ! Standard value = false = missing solutions are set to zero.
  LOGICAL, OPTIONAL :: brepeatLast

  ! OPTIONAL: Directory name where to save temporary data. If not present,
  ! the current directory './' is the default.
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sdirectory
  
!</input>

!<inputoutput>
  ! Space-time vector where store data to.
  ! The vector is created from the scratch.
  TYPE(t_spaceTimeVector), INTENT(OUT) :: rx
!</inputoutput>

!</subroutine>

    ! Local variables
    TYPE(t_vectorScalar) :: rvector
    CHARACTER(SYS_STRLEN) :: sfile,sarray
    INTEGER :: i,ilast
    LOGICAL :: bexists,brepeat
    
    brepeat = .FALSE.
    IF (PRESENT(brepeatLast)) brepeat = brepeatLast
    
    ilast = 0

    ! Loop over the files
    DO i=istart,iend
    
      ! Form the filename
      WRITE(sfile,sfilename) i*idelta
      
      ! Is the file there?
      INQUIRE(file=trim(sfile), exist=bexists)
      
      IF (bexists) THEN
        ! Remember this solution as the last available one
        ilast = i
      
        ! Read the file into rvector. The first read command creates rvector
        ! in the correct size.
        CALL vecio_readVectorHR (rvector, sarray, .FALSE.,&
          0, sfile, bformatted)

        IF (i .EQ. istart) THEN
          ! At the first file, create a space-time vector holding the data.
          CALL sptivec_initVector (rx,rvector%NEQ,iend-istart+1,sdirectory)
        END IF         
        
        ! Save the data
        CALL storage_copy (rvector%h_Ddata,rx%p_IdataHandleList(i))
      ELSE
        IF (i .EQ. istart) THEN
          ! At the first file, create a space-time vector holding the data.
          CALL sptivec_initVector (rx,rvector%NEQ,iend-istart+1,sdirectory)
        END IF         

        IF (brepeat) THEN
          CALL output_line ('Warning: Unable to load file "'//TRIM(sfile) &
              //'". Repeating last solution!', &
              ssubroutine='sptivec_loadFromFileSequence')
        
          ! Copy the data from the last known solution to the current one.
          CALL storage_copy (&
              rx%p_IdataHandleList(ilast),&
              rx%p_IdataHandleList(i))
        ELSE
          CALL output_line ('Warning: Unable to load file "'//TRIM(sfile) &
              //'". Assuming zero!', ssubroutine='sptivec_loadFromFileSequence')
        
          ! Clear that array. Zero solution.
          CALL storage_clear (rx%p_IdataHandleList(i))
        END IF
      END IF
    
    END DO
    
    ! The vector is scaled by 1.0.
    rx%p_Dscale(:) = 1.0_DP
    
    ! Remove the temp vector
    CALL lsyssc_releaseVector (rvector)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_saveToFileSequence (rx,sfilename,bformatted,&
      rtempVector,sformat)

!<description>
  ! This routine saves a space-time vector to a sequence of files on the
  ! hard disc. sfilename is a directory/file name pattern in the format of
  ! a FORMAT statement that forms the filename; this pattern must contain
  ! exactly one integer format specifier, which is replaced by the
  ! file number in this routine (e.g. ' (''vector.txt.'',I5.5) ' will
  ! load a file sequence 'vector.txt.00000','vector.txt.00001',
  ! 'vector.txt.00002', ...).
  !
  ! The destination files are written out using vecio_writeBlockVectorHR 
  ! or vecio_writeVectorHR. Old files are overwritten.
!</description>

!<input>
  ! Space-time vector to be written to disc.
  TYPE(t_spaceTimeVector), INTENT(IN) :: rx

  ! Filename pattern + path where to form a filename from.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename

  ! Whether to read formatted or unformatted data from disc.
  LOGICAL, INTENT(IN) :: bformatted
  
  ! OPTIONAL: Format string that is used for exporting data to files.
  ! E.g. '(E20.10)'.
  CHARACTER(LEN=SYS_STRLEN), INTENT(IN), OPTIONAL :: sformat
!</input>

!<inputoutput>
  ! Temporary vector. This vector must prescribe the block structure of the
  ! subvectors in the space-time vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector
!</inputoutput>

!</subroutine>

    ! Local variables
    CHARACTER(SYS_STRLEN) :: sfile
    INTEGER :: i
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx

    CALL lsysbl_getbase_double (rtempVector,p_Dx)

    ! Loop over the files
    DO i=0,rx%ntimesteps
    
      ! Get the data from the space-time vector.
      CALL sptivec_getTimestepData (rx, i, rtempVector)
      
      ! Form the filename
      WRITE(sfile,sfilename) i
      
      ! Save that to disc.
      IF (.NOT. bformatted) THEN
        CALL vecio_writeBlockVectorHR (rtempVector, 'vector', .TRUE.,&
                                      0, sfile)
      ELSE IF (.NOT. PRESENT(sformat)) THEN
          CALL vecio_writeBlockVectorHR (rtempVector, 'vector', .TRUE.,&
                                        0, sfile, '(E24.16)')
      ELSE
        CALL vecio_writeBlockVectorHR (rtempVector, 'vector', .TRUE.,&
                                       0, sfile, sformat)
      END IF
    
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_setConstant (rx,dvalue)

!<description>
  ! Sets the whole space-time vector to a constant value.
!</description>

!<input>
  ! Value which should be written to the whole space-time vector.
  REAL(DP), INTENT(IN) :: dvalue
!</input>

!<inputoutput>
  ! Space-time vector to modify
  TYPE(t_spaceTimeVector), INTENT(INOUT) :: rx
!</inputoutput>

!</subroutine>

    ! Local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER :: i

    ! Loop over the files
    DO i=0,rx%ntimesteps
    
      ! Get the data and set to a defined value.
      CALL storage_getbase_double (rx%p_IdataHandleList(i),p_Ddata)
      p_Ddata(:) = dvalue
      rx%p_Dscale(i) = 1.0_DP

    END DO
    
  END SUBROUTINE

END MODULE
