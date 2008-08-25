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
!# Making use of the externalstorage.f90 library, vectors can be held on
!# external memory storage devices.
!#
!# The module provides the following subroutines:
!#
!# 1.) sptivec_initVector = 
!#     sptivec_initVectorPlain, 
!#     sptivec_initVectorDirect,
!#     sptivec_initVectorDiscr
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
!# 17.) sptivec_scalarProduct
!#     -> Calculate the scalar product of two vectors
!#
!# 18.) sptivec_scalarProductWeighted
!#     -> Calculate the weighted scalar product of two vectors
!#
!# </purpose>
!##############################################################################

MODULE spacetimevectors

  USE fsystem
  USE genoutput
  USE externalstorage
  USE spatialdiscretisation
  USE linearsystemscalar
  USE linearsystemblock
  USE collection
  USE vectorio
  
  USE timediscretisation

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
  
    ! This flag defines a scaling factor for each substep. The scaling is applied
    ! when a subvector is read and is reset to 1.0 if a subvector is saved.
    ! The whole vector is zero if it's new and if it's cleared by clearVector.
    REAL(DP), DIMENSION(:), POINTER :: p_Dscale => NULL()
    
    ! Pointer to a time discretisation structure that defines the
    ! discretisation in time.
    TYPE(t_timeDiscretisation), POINTER :: p_rtimeDiscretisation => NULL()
    
    ! Pointer to the underlying block discretisation that defines the
    ! spatial discretisation.
    TYPE(t_blockDiscretisation), POINTER :: p_rblockDiscretisation => NULL()
    
    ! Whether to use the test space of the spatial discretisation for the
    ! shape of the spatial vectors or nor.
    LOGICAL :: btestfctSpace = .FALSE.
    
    ! Number of equations in each subvector of the 'global time-step vector'.
    INTEGER(PREC_VECIDX) :: NEQ = 0
    
    ! Number of subvectors saved in p_IdataHandleList.
    INTEGER :: NEQtime = 0
    
    ! A list of handles (dimension 1:NEQtime) to double-precision arrays 
    ! which save the data of the ntimesteps+1 data subvectors. 
    ! A value ST_NOHANDLE indicates that there is currently no data stored
    ! at that timestep.
    INTEGER, DIMENSION(:), POINTER :: p_IdataHandleList => NULL()
    
  END TYPE

!</typeblock>

!</types>

  INTERFACE sptivec_initVector
    MODULE PROCEDURE sptivec_initVectorPlain
    MODULE PROCEDURE sptivec_initVectorDirect 
    MODULE PROCEDURE sptivec_initVectorDiscr
  END INTERFACE
    
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_initVectorPlain (rspaceTimeVector,NEQ,NEQtime)

!<description>
  ! Initialises a space time vector. NEQ defines the size of each spatial
  ! subvector. ntimesteps defines the number of timesteps to maintain.
  ! The subvectors have no special block structure, they are just plain
  ! data.
!</desctiprion>

!<input>
  ! Number of equations in the vectors
  INTEGER(PREC_VECIDX), INTENT(IN) :: NEQ
  
  ! Number of DOF's in time to maintain.
  ! The number of subvectors that is reserved is therefore ntimesteps+1!
  INTEGER, INTENT(IN) :: NEQtime
!</input>

!<output>
  ! Space-time vector structure to be initialised.
  TYPE(t_spacetimeVector), INTENT(OUT) :: rspaceTimeVector
!</output>

!</subroutine>

    INTEGER :: i

    ! Initialise the data.
    rspaceTimeVector%NEQtime = NEQtime
    ALLOCATE(rspaceTimeVector%p_IdataHandleList(1:NEQtime))
    ALLOCATE(rspaceTimeVector%p_Dscale(1:NEQtime))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    rspaceTimeVector%p_Dscale(:) = 0.0_DP
    
    rspaceTimeVector%NEQ = NEQ
    
    ! Allocate memory for every subvector
    DO i=1,NEQtime
      CALL exstor_new ('sptivec_initVector', 'stvec_'//TRIM(sys_siL(i,10)), &
        NEQ, ST_DOUBLE, rspaceTimeVector%p_IdataHandleList(i), ST_NEWBLOCK_ZERO)
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_initVectorDirect (rspaceTimeVector,NEQtime,rblockDiscr)

!<description>
  ! Initialises a space time vector. rblockDiscr is a block discretisation 
  ! structure that defines the basic shape of the data vectors in all time 
  ! steps that are to be maintained. ntimesteps defines the number of 
  ! timesteps to maintain.
!</desctiprion>

!<input>
  ! Number of DOF's in time to maintain.
  ! The number of subvectors that is reserved is therefore ntimesteps+1!
  INTEGER, INTENT(IN) :: NEQtime

  ! Block discretisation structure of the spatial discretisation.
  ! A pointer to this structure is saved in the space time vector.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: rblockDiscr
  
!</input>

!<output>
  ! Space-time vector structure to be initialised.
  TYPE(t_spacetimeVector), INTENT(OUT) :: rspaceTimeVector
!</output>

!</subroutine>

    INTEGER :: i

    ! Initialise the data.
    rspaceTimeVector%NEQtime = NEQtime
    ALLOCATE(rspaceTimeVector%p_IdataHandleList(1:NEQtime))
    ALLOCATE(rspaceTimeVector%p_Dscale(1:NEQtime))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    rspaceTimeVector%p_Dscale(:) = 0.0_DP
    
    ! Get NEQ and save a pointer to the spatial discretisation structure
    ! to the vector.
    rspaceTimeVector%p_rblockDiscretisation => rblockDiscr
    rspaceTimeVector%NEQ = dof_igetNDofGlobBlock(rblockDiscr)
    
    ! Allocate memory for every subvector
    DO i=1,NEQtime
      CALL exstor_new ('sptivec_initVector', 'stvec_'//TRIM(sys_siL(i,10)), &
        rspaceTimeVector%NEQ, ST_DOUBLE, rspaceTimeVector%p_IdataHandleList(i), &
        ST_NEWBLOCK_ZERO)
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE sptivec_initVectorDiscr (rspaceTimeVector,rtimeDiscr,rblockDiscr)

!<description>
  ! Initialises a space time vector according to a time discretisation
  ! structure. rblockDiscr is a block discretisation 
  ! structure that defines the basic shape of the data vectors in all time 
  ! steps that are to be maintained. 
!</desctiprion>

!<input>

  ! Time discretisation structure that defines the discrtetisation
  ! in time. A pointer to this is saved to rspaceTimeVector.
  TYPE(t_timeDiscretisation), INTENT(IN), TARGET :: rtimeDiscr

  ! Block discretisation structure of the spatial discretisation.
  ! A pointer to this structure is saved in the space time vector.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: rblockDiscr
  
!</input>

!<output>
  ! Space-time vector structure to be initialised.
  TYPE(t_spacetimeVector), INTENT(OUT) :: rspaceTimeVector
!</output>

!</subroutine>

    INTEGER :: i,NEQtime

    ! Initialise the data.

    ! Get the number of DOF's in time.    
    NEQtime = tdiscr_igetNDofGlob(rtimediscr)
    
    rspaceTimeVector%NEQtime = NEQtime
    rspaceTimeVector%p_rtimeDiscretisation => rtimeDiscr 
    
    ALLOCATE(rspaceTimeVector%p_IdataHandleList(1:NEQtime))
    ALLOCATE(rspaceTimeVector%p_Dscale(1:NEQtime))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    rspaceTimeVector%p_Dscale(:) = 0.0_DP
    
    ! Get NEQ and save a pointer to the spatial discretisation structure
    ! to the vector.
    rspaceTimeVector%p_rblockDiscretisation => rblockDiscr
    rspaceTimeVector%NEQ = dof_igetNDofGlobBlock(rblockDiscr)
    
    ! Allocate memory for every subvector
    DO i=1,NEQtime
      CALL exstor_new ('sptivec_initVector', 'stvec_'//TRIM(sys_siL(i,10)), &
        rspaceTimeVector%NEQ, ST_DOUBLE, rspaceTimeVector%p_IdataHandleList(i), &
        ST_NEWBLOCK_ZERO)
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
      DO i=1,rspaceTimeVector%NEQtime
        IF (rspaceTimeVector%p_IdataHandleList(i) .NE. ST_NOHANDLE) THEN
        
          CALL exstor_free (rspaceTimeVector%p_IdataHandleList(i))
          
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
  ! Number of the subvector that corresponds to rvector. >= 1, <= NEQtime from
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
    REAL(DP), DIMENSION(:), POINTER :: p_Dsource

    ! Make sure we can store the timestep data.
    IF ((isubvector .LT. 1) .OR. (isubvector .GT. rspaceTimeVector%NEQtime)) THEN
      CALL output_line('Invalid timestep number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_setTimestepData')
      CALL sys_halt()
    END IF
    
    IF (rvector%NEQ .NE. rspaceTimeVector%NEQ) THEN
      CALL output_line('Vector size invalid!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_setTimestepData')
      CALL sys_halt()
    END IF

    ! Save the vector data. If necessary, new memory is allocated -- as the
    ! default value of the handle in the handle list is ST_NOHANDLE.
    !
    ! Don't use storage_copy, as this might give errors in case the array
    ! behind the handle is longer than the vector!
    CALL lsysbl_getbase_double (rvector,p_Dsource)
    !CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(isubvector),p_Ddest)
    !CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
    CALL exstor_setdata_double (rspaceTimeVector%p_IdataHandleList(isubvector),p_Dsource)

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
  ! Number of the subvector that corresponds to rvector. >= 1, <= NEQtime from
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
    REAL(DP), DIMENSION(:), POINTER :: p_Ddest

    ! Make sure we can store the timestep data.
    IF ((isubvector .LT. 1) .OR. (isubvector .GT. rspaceTimeVector%NEQtime)) THEN
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

    ! Get the vector data. 
    !
    ! Don't use storage_copy, as this might give errors in case the array
    ! behind the handle is longer than the vector!
    !CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(isubvector),&
    !    p_Dsource)
    CALL lsysbl_getbase_double (rvector,p_Ddest)
    !CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
    CALL exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(isubvector),p_Ddest)
    
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
  ! stamp in the range 0.0 .. 1.0, where 0.0 corresponds to the 1st subvector
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
    REAL(DP), DIMENSION(:), POINTER :: p_Ddest
    INTEGER :: itimestep1,itimestep2,itimestep3
    REAL(DP) :: dreltime,dabstime
    INTEGER :: i
    REAL(DP) :: dscal1,dscal2,dscal3
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_Dsource

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
    ! Rescale dtimestamp to the interval [1..NEQtime].
    dabstime = dtimestamp*REAL(rspaceTimeVector%NEQtime-1,DP)+1.0_DP
    itimestep2 = INT(dabstime + 0.5_DP)
    
    IF (dabstime .EQ. REAL(itimestep2,DP)) THEN
      ! Nice coincidence, we have exactly timestep itimestep2. Ok, then we 
      ! can call the routine to get that timestep; this saves us some
      ! time as the interpolation can be omitted.
      CALL sptivec_getTimestepData (rspaceTimeVector, itimestep2, rvector)
      RETURN
    END IF
    
    IF (rspaceTimeVector%NEQtime .EQ. 2) THEN
      ! Special case: only one timestep!
      itimestep1 = 0
      itimestep2 = 0
      itimestep3 = 1
    ELSE
      ! Is this the first or the last timestep?
      IF (itimestep2 .EQ. 1) THEN
        ! First timestep. Interpolate between timesteps 0,1 and 2, evaluate 
        ! near timestep 0.
        itimestep1 = 0
        itimestep2 = 1
        itimestep3 = 2
      ELSE IF (itimestep2 .EQ. rspaceTimeVector%NEQtime) THEN
        ! Last timestep. Interpolate between timesteps n-2,n-1 and n, evaluate 
        ! near timestep n.
        itimestep1 = rspaceTimeVector%NEQtime-2
        itimestep2 = rspaceTimeVector%NEQtime-1
        itimestep3 = rspaceTimeVector%NEQtime
      ELSE
        ! Somewhere in the inner. Get the number of the previous and next timestep
        itimestep1 = itimestep2-1
        itimestep3 = itimestep2+1
      END IF
    END IF

    ! Calculate the 'relative' time in the interval [-1,1], where -1 corresponds
    ! to timestep itimestep1, 0 to itimestep2 and +1 to itimestep3.
    ! This will be used to evaluate the quadratic polynomial.
    dreltime = dabstime-1.0_DP-REAL(itimestep2,DP)
    
    ! Get the vector data of the three timesteps
    ALLOCATE(p_Dsource(rspaceTimeVector%NEQ,3))
    
    CALL exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(itimestep1),&
        p_Dsource(:,1))
    CALL exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(itimestep2),&
        p_Dsource(:,2))
    CALL exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(itimestep3),&
        p_Dsource(:,3))

    CALL lsysbl_getbase_double (rvector,p_Ddest)

    ! Calculate the quadratic interpolation of the three arrays.
    dscal1 = rspaceTimeVector%p_Dscale(itimestep1)
    dscal2 = rspaceTimeVector%p_Dscale(itimestep2)
    dscal3 = rspaceTimeVector%p_Dscale(itimestep3)
    IF (rspaceTimeVector%NEQtime .EQ. 2) THEN
      ! Special case: only 1 timestep. Linear interpolation. dreltime is in [0..1]!
      DO i=1,SIZE(p_Ddest)
        p_Ddest(i) = (1.0_DP-dreltime) * dscal2*p_Dsource(i,2) + &
                     dreltime*dscal3*p_Dsource(i,3)
      END DO
    ELSE
      ! Quadratic interpolation
      DO i=1,SIZE(p_Ddest)
        CALL mprim_quadraticInterpolation (dreltime,&
            dscal1*p_Dsource(i,1),dscal2*p_Dsource(i,2),dscal3*p_Dsource(i,3),p_Ddest(i))
      END DO
    END IF
    
    DEALLOCATE(p_Dsource)
    
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
      CALL output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_vectorLinearComb')
      CALL sys_halt()
    END IF

    IF (rx%NEQtime .NE. ry%NEQtime) THEN
      CALL output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_vectorLinearComb')
      CALL sys_halt()
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
    DO i=1,rx%NEQtime
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
    
    IF ((ry%NEQ .EQ. 0) .AND. (ry%NEQtime .EQ. 0)) THEN
      ! Destination vector does not exist. Create it.
      CALL sptivec_initVector (ry,rx%NEQtime,rx%p_rblockDiscretisation)
    END IF
    
    IF (rx%NEQ .NE. ry%NEQ) THEN
      CALL output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_copyVector')
      CALL sys_halt()
    END IF

    IF (rx%NEQtime .NE. ry%NEQtime) THEN
      CALL output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_copyVector')
      CALL sys_halt()
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
    DO i=1,rx%NEQtime
      
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
      CALL output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      CALL sys_halt()
    END IF

    IF (rx%NEQtime .NE. ry%NEQtime) THEN
      CALL output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      CALL sys_halt()
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
    DO i=1,rx%NEQtime
      
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
  
  REAL(DP) FUNCTION sptivec_scalarProductWeighted (rx, ry, Dweights)
  
!<description>
  ! Calculates a weighted scalar product of two block vectors.
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>

!<input>
  ! First source vector
  TYPE(t_spacetimeVector), INTENT(IN)   :: rx

  ! Second source vector
  TYPE(t_spacetimeVector), INTENT(IN)   :: ry
  
  ! An array with weights for each component of a solution
  ! vector in a timestep. Dweights(i) is multiplied to
  ! the i'th solution component of each timestep.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dweights
!</input>

!<result>
  ! The scalar product (rx,ry) of the two block vectors.
!</result>

!</function>

    INTEGER :: i,irow
    INTEGER(PREC_VECIDX), DIMENSION(1) :: Isize
    REAL(DP) :: dres,a
    TYPE(t_vectorBlock) :: rxBlock,ryBlock
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dy
    
    IF (rx%NEQ .NE. ry%NEQ) THEN
      CALL output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      CALL sys_halt()
    END IF

    IF (rx%NEQtime .NE. ry%NEQtime) THEN
      CALL output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      CALL sys_halt()
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
    DO i=1,rx%NEQtime
      
      IF ((rx%p_Dscale(i) .NE. 0.0_DP) .AND. (ry%p_Dscale(i) .NE. 0.0_DP)) THEN
        CALL sptivec_getTimestepData (rx, i, rxBlock)
        CALL sptivec_getTimestepData (ry, i, ryBlock)
        
        ! Calculate a weighted scalar product using the equation weights.
        a = 0.0_DP
        DO irow = 1,rxBlock%nblocks
          a = a + Dweights(irow) * &
              lsyssc_scalarProduct(rxBlock%RvectorBlock(irow),&
                                   ryBlock%RvectorBlock(irow))
        END DO
        
        dres = dres + a
      END IF

    END DO

    ! Release temp memory    
    CALL lsysbl_releaseVector (ryBlock)
    CALL lsysbl_releaseVector (rxBlock)
      
    sptivec_scalarProductWeighted = dres
      
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
    DO i=1,rx%NEQtime
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
      sptivec_vectorNorm = dnorm / (rx%NEQtime)
    CASE (LINALG_NORML2)
      sptivec_vectorNorm = SQRT(dnorm / (rx%NEQtime))
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
    ! rx%p_Dscale(:) = 0.0_DP
    CALL lalg_clearVectorDble (rx%p_Dscale(:))

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
    ! rx%p_Dscale(:) = rx%p_Dscale(:) * dscale
    CALL lalg_scaleVectorDble(rx%p_Dscale(:),dscale)

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
    ALLOCATE(Isize(rxsuper%NEQtime))
    Isize(:) = rxsuper%NEQ
    CALL lsysbl_createVecBlockDirect (rx,Isize,.FALSE.)

    ! Create a 1-block temp vector for the data    
    CALL lsysbl_createVecBlockDirect (rvectorTmp,Isize(1:1),.FALSE.)
    CALL lsysbl_getbase_double (rvectorTmp,p_Ddata1)
    
    ! Load the subvectors and write them to the global vector.
    DO i=1,rxsuper%NEQtime
      CALL sptivec_getTimestepData (rxsuper, i, rvectorTmp)
      
      CALL lsyssc_getbase_double (rx%RvectorBlock(i),p_Ddata2)
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
    DO i=1,rxsuper%NEQtime
      CALL lsyssc_getbase_double (rx%RvectorBlock(i),p_Ddata1)
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
    CALL collct_setvalue_int (rcollection,TRIM(sname)//'_NEQ',&
        rx%NEQ,.TRUE.,ilevel,ssection)
  
    CALL collct_setvalue_int (rcollection,TRIM(sname)//'_NTST',&
        rx%NEQtime,.TRUE.,ilevel,ssection)

    IF (rx%NEQtime .NE. 0) THEN
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
    rx%NEQ = collct_getvalue_int (rcollection,TRIM(sname)//'_NEQ',&
        ilevel,ssection)
  
    rx%NEQtime = collct_getvalue_int (rcollection,TRIM(sname)//'_NTST',&
        ilevel,ssection)
        
    IF (rx%NEQtime .NE. 0) THEN
      ! For the handle list, we need to allocate some memory...
      ALLOCATE(rx%p_IdataHandleList(rx%NEQtime))

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
      bformatted,brepeatLast,rblockDiscretisation)

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
  
  ! OPTIONAL: A block discretisation structure that defines the shape of the spatial
  ! vectors. If not specified, the vectors will be read in as pure data vectors
  ! without a discretisation attached.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET, OPTIONAL :: rblockDiscretisation

!</input>

!<inputoutput>
  ! Space-time vector where store data to.
  ! If empty, the vector is created from the scratch.
  ! If the vector is not empty, it must be large enough to hold all data.
  ! The content is ovrwritten.
  TYPE(t_spaceTimeVector), INTENT(INOUT) :: rx
!</inputoutput>

!</subroutine>

    ! Local variables
    TYPE(t_vectorBlock) :: rvector
    TYPE(t_vectorScalar) :: rvectorScalar
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
      
        IF (.NOT. PRESENT(rblockDiscretisation)) THEN
          ! Read the file into rvector. The first read command creates rvector
          ! in the correct size.
          CALL vecio_readVectorHR (rvectorScalar, sarray, .FALSE.,&
            0, sfile, bformatted)

          IF ((i .EQ. istart) .AND. (rx%NEQtime .EQ. 0)) THEN
            ! At the first file, create a space-time vector holding the data.
            CALL sptivec_initVectorPlain (rx,rvectorScalar%NEQ,iend-istart+1)
          END IF         

          ! Save the data
          CALL exstor_setdata_storage (rx%p_IdataHandleList(1+i),rvectorScalar%h_Ddata)
    
        ELSE
  
          CALL vecio_readBlockVectorHR (rvector, sarray, .FALSE.,&
            0, sfile, bformatted)

          IF (i .EQ. istart) THEN
            ! At the first file, create a space-time vector holding the data.
            CALL sptivec_initVector (rx,1+iend-istart,rblockDiscretisation)
          END IF         
          ! Save the data
          CALL exstor_setdata_storage (rx%p_IdataHandleList(1+i),rvector%h_Ddata)
          
        END IF
        
      ELSE
        IF (i .EQ. istart) THEN
          ! The first file must exist!
          CALL output_line ('The first file must exist!', &
              ssubroutine='sptivec_loadFromFileSequence')
        END IF         

        IF (brepeat) THEN
          CALL output_line ('Warning: Unable to load file "'//TRIM(sfile) &
              //'". Repeating last solution!', &
              ssubroutine='sptivec_loadFromFileSequence')
        
          ! Copy the data from the last known solution to the current one.
          CALL exstor_copy (&
              rx%p_IdataHandleList(1+ilast),&
              rx%p_IdataHandleList(1+i))
        ELSE
          CALL output_line ('Warning: Unable to load file "'//TRIM(sfile) &
              //'". Assuming zero!', ssubroutine='sptivec_loadFromFileSequence')
        
          ! Clear that array. Zero solution.
          CALL exstor_clear (rx%p_IdataHandleList(1+i))
        END IF
      END IF
    
    END DO
    
    ! The vector is scaled by 1.0.
    ! rx%p_Dscale(:) = 1.0_DP
    CALL lalg_setVectorDble(rx%p_Dscale(:),1.0_DP)
    
    ! Remove the temp vector
    IF (PRESENT(rblockDiscretisation)) THEN
      CALL lsysbl_releaseVector (rvector)
    ELSE
      CALL lsyssc_releaseVector (rvectorScalar)
    END IF

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
  ! subvectors in the space-time vector. If not specified, the data is written
  ! out without a block structure.
  TYPE(t_vectorBlock), INTENT(INOUT), TARGET, OPTIONAL :: rtempVector
!</inputoutput>

!</subroutine>

    ! Local variables
    CHARACTER(SYS_STRLEN) :: sfile
    INTEGER :: i
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx
    TYPE(t_vectorBlock), POINTER :: p_rx
    INTEGER(PREC_VECIDX), DIMENSION(1) :: Isize
    
    IF (PRESENT(rtempVector)) THEN
      p_rx => rtempVector
    ELSE
      ! Create a 1-block temp vector for the data    
      ALLOCATE(p_rx)
      Isize(1) = rx%NEQ
      CALL lsysbl_createVecBlockDirect (p_rx,Isize(1:1),.FALSE.)
    END IF

    ! DEBUG!!!
    CALL lsysbl_getbase_double (p_rx,p_Dx)

    ! Loop over the files
    DO i=0,rx%NEQtime-1
    
      ! Get the data from the space-time vector.
      CALL sptivec_getTimestepData (rx, 1+i, p_rx)
      
      ! Form the filename
      WRITE(sfile,sfilename) i
      
      ! Save that to disc.
      IF (.NOT. bformatted) THEN
        CALL vecio_writeBlockVectorHR (p_rx, 'vector', .TRUE.,&
                                      0, sfile)
      ELSE IF (.NOT. PRESENT(sformat)) THEN
          CALL vecio_writeBlockVectorHR (p_rx, 'vector', .TRUE.,&
                                        0, sfile, '(E24.16)')
      ELSE
        CALL vecio_writeBlockVectorHR (p_rx, 'vector', .TRUE.,&
                                       0, sfile, sformat)
      END IF
    
    END DO
    
    IF (.NOT. PRESENT(rtempVector)) THEN
      CALL lsysbl_releaseVector (p_rx)
      DEALLOCATE(p_rx)
    END IF
    
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
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_Ddata
    INTEGER :: i
    
    ! Allocate memory for intermediate values
    ALLOCATE(p_Ddata(rx%NEQ))

    ! Loop over the files
    DO i=1,rx%NEQtime
    
      ! Get the data and set to a defined value.
      CALL exstor_getdata_double (rx%p_IdataHandleList(i),p_Ddata)
      p_Ddata(:) = dvalue
      CALL exstor_setdata_double (rx%p_IdataHandleList(i),p_Ddata)

      rx%p_Dscale(i) = 1.0_DP

    END DO
    
    DEALLOCATE(p_Ddata)
    
  END SUBROUTINE

END MODULE
