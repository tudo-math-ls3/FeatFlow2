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

module spacetimevectors

  use fsystem
  use genoutput
  use externalstorage
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  
  use timediscretisation

  implicit none

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
  type t_spacetimeVector
  
    ! Whether this vector shares its data with another vector.
    logical :: bisCopy = .false.
  
    ! This flag defines a scaling factor for each substep. The scaling is applied
    ! when a subvector is read and is reset to 1.0 if a subvector is saved.
    ! The whole vector is zero if it's new and if it's cleared by clearVector.
    real(DP), dimension(:), pointer :: p_Dscale => null()
    
    ! Pointer to a time discretisation structure that defines the
    ! discretisation in time.
    type(t_timeDiscretisation), pointer :: p_rtimeDiscretisation => null()
    
    ! Pointer to the underlying block discretisation that defines the
    ! spatial discretisation.
    type(t_blockDiscretisation), pointer :: p_rblockDiscretisation => null()
    
    ! Whether to use the test space of the spatial discretisation for the
    ! shape of the spatial vectors or nor.
    logical :: btestfctSpace = .false.
    
    ! Number of equations in each subvector of the 'global time-step vector'.
    integer(PREC_VECIDX) :: NEQ = 0
    
    ! Number of subvectors saved in p_IdataHandleList.
    integer :: NEQtime = 0
    
    ! A list of handles (dimension 1:NEQtime) to double-precision arrays 
    ! which save the data of the ntimesteps+1 data subvectors. 
    ! A value ST_NOHANDLE indicates that there is currently no data stored
    ! at that timestep.
    integer, dimension(:), pointer :: p_IdataHandleList => null()
    
  end type

!</typeblock>

!</types>

  interface sptivec_initVector
    module procedure sptivec_initVectorPlain
    module procedure sptivec_initVectorDirect 
    module procedure sptivec_initVectorDiscr
  end interface
    
contains

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_initVectorPlain (rspaceTimeVector,NEQ,NEQtime)

!<description>
  ! Initialises a space time vector. NEQ defines the size of each spatial
  ! subvector. ntimesteps defines the number of timesteps to maintain.
  ! The subvectors have no special block structure, they are just plain
  ! data.
!</desctiprion>

!<input>
  ! Number of equations in the vectors
  integer(PREC_VECIDX), intent(IN) :: NEQ
  
  ! Number of DOF's in time to maintain.
  ! The number of subvectors that is reserved is therefore ntimesteps+1!
  integer, intent(IN) :: NEQtime
!</input>

!<output>
  ! Space-time vector structure to be initialised.
  type(t_spacetimeVector), intent(OUT) :: rspaceTimeVector
!</output>

!</subroutine>

    integer :: i

    ! Initialise the data.
    rspaceTimeVector%NEQtime = NEQtime
    allocate(rspaceTimeVector%p_IdataHandleList(1:NEQtime))
    allocate(rspaceTimeVector%p_Dscale(1:NEQtime))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    rspaceTimeVector%p_Dscale(:) = 0.0_DP
    
    rspaceTimeVector%NEQ = NEQ
    
    ! Allocate memory for every subvector
    do i=1,NEQtime
      call exstor_new ('sptivec_initVector', 'stvec_'//trim(sys_siL(i,10)), &
        NEQ, ST_DOUBLE, rspaceTimeVector%p_IdataHandleList(i), ST_NEWBLOCK_ZERO)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_initVectorDirect (rspaceTimeVector,NEQtime,rblockDiscr)

!<description>
  ! Initialises a space time vector. rblockDiscr is a block discretisation 
  ! structure that defines the basic shape of the data vectors in all time 
  ! steps that are to be maintained. ntimesteps defines the number of 
  ! timesteps to maintain.
!</desctiprion>

!<input>
  ! Number of DOF's in time to maintain.
  ! The number of subvectors that is reserved is therefore ntimesteps+1!
  integer, intent(IN) :: NEQtime

  ! Block discretisation structure of the spatial discretisation.
  ! A pointer to this structure is saved in the space time vector.
  type(t_blockDiscretisation), intent(IN), target :: rblockDiscr
  
!</input>

!<output>
  ! Space-time vector structure to be initialised.
  type(t_spacetimeVector), intent(OUT) :: rspaceTimeVector
!</output>

!</subroutine>

    integer :: i

    ! Initialise the data.
    rspaceTimeVector%NEQtime = NEQtime
    allocate(rspaceTimeVector%p_IdataHandleList(1:NEQtime))
    allocate(rspaceTimeVector%p_Dscale(1:NEQtime))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    rspaceTimeVector%p_Dscale(:) = 0.0_DP
    
    ! Get NEQ and save a pointer to the spatial discretisation structure
    ! to the vector.
    rspaceTimeVector%p_rblockDiscretisation => rblockDiscr
    rspaceTimeVector%NEQ = dof_igetNDofGlobBlock(rblockDiscr)
    
    ! Allocate memory for every subvector
    do i=1,NEQtime
      call exstor_new ('sptivec_initVector', 'stvec_'//trim(sys_siL(i,10)), &
        rspaceTimeVector%NEQ, ST_DOUBLE, rspaceTimeVector%p_IdataHandleList(i), &
        ST_NEWBLOCK_ZERO)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_initVectorDiscr (rspaceTimeVector,rtimeDiscr,rblockDiscr)

!<description>
  ! Initialises a space time vector according to a time discretisation
  ! structure. rblockDiscr is a block discretisation 
  ! structure that defines the basic shape of the data vectors in all time 
  ! steps that are to be maintained. 
!</desctiprion>

!<input>

  ! Time discretisation structure that defines the discrtetisation
  ! in time. A pointer to this is saved to rspaceTimeVector.
  type(t_timeDiscretisation), intent(IN), target :: rtimeDiscr

  ! Block discretisation structure of the spatial discretisation.
  ! A pointer to this structure is saved in the space time vector.
  type(t_blockDiscretisation), intent(IN), target :: rblockDiscr
  
!</input>

!<output>
  ! Space-time vector structure to be initialised.
  type(t_spacetimeVector), intent(OUT) :: rspaceTimeVector
!</output>

!</subroutine>

    integer :: i,NEQtime

    ! Initialise the data.

    ! Get the number of DOF's in time.    
    NEQtime = tdiscr_igetNDofGlob(rtimediscr)
    
    rspaceTimeVector%NEQtime = NEQtime
    rspaceTimeVector%p_rtimeDiscretisation => rtimeDiscr 
    
    allocate(rspaceTimeVector%p_IdataHandleList(1:NEQtime))
    allocate(rspaceTimeVector%p_Dscale(1:NEQtime))
    rspaceTimeVector%p_IdataHandleList(:) = ST_NOHANDLE
    rspaceTimeVector%p_Dscale(:) = 0.0_DP
    
    ! Get NEQ and save a pointer to the spatial discretisation structure
    ! to the vector.
    rspaceTimeVector%p_rblockDiscretisation => rblockDiscr
    rspaceTimeVector%NEQ = dof_igetNDofGlobBlock(rblockDiscr)
    
    ! Allocate memory for every subvector
    do i=1,NEQtime
      call exstor_new ('sptivec_initVector', 'stvec_'//trim(sys_siL(i,10)), &
        rspaceTimeVector%NEQ, ST_DOUBLE, rspaceTimeVector%p_IdataHandleList(i), &
        ST_NEWBLOCK_ZERO)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_releaseVector (rspaceTimeVector)

!<description>
  ! Releases a space-time vector. All allocated memory is released. Temporary data
  ! on disc is deleted.
!</desctiprion>

!<inputoutput>
  ! Space-time vector structure to be initialised.
  type(t_spacetimeVector), intent(INOUT) :: rspaceTimeVector
!</inputoutput>

!</subroutine>

    ! local variables -- initialised by Fortran default initialisation!
    type(t_spacetimeVector) :: rspaceTimeVectorTempl
    integer :: i

    if (.not. associated(rspaceTimeVector%p_IdataHandleList)) then
      call output_line('Warning: Releasing unused vector!',&
          ssubroutine='sptivec_releaseVector')
      return
    end if

    ! Deallocate data -- if the data is not shared with another vector
    if (.not. rspaceTimeVector%bisCopy) then
      do i=1,rspaceTimeVector%NEQtime
        if (rspaceTimeVector%p_IdataHandleList(i) .ne. ST_NOHANDLE) then
        
          call exstor_free (rspaceTimeVector%p_IdataHandleList(i))
          
        end if
      end do
    end if

    deallocate(rspaceTimeVector%p_IdataHandleList)
    deallocate(rspaceTimeVector%p_Dscale)

    ! Initialise with default values.
    rspaceTimeVector = rspaceTimeVectorTempl

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_setTimestepData (rspaceTimeVector, isubvector, rvector)

!<description>
  ! Stores the data of rvector at timestep itimestep in the rspaceTimeVector.
  ! structure.
!</desctiprion>

!<input>
  ! Number of the subvector that corresponds to rvector. >= 1, <= NEQtime from
  ! the initialisation of rspaceTimeVector.
  integer, intent(IN) :: isubvector
  
  ! Vector with data that should be associated to timestep itimestep.
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!<inputoutput>
  ! Space-time vector structure where to save the data.
  type(t_spacetimeVector), intent(INOUT) :: rspaceTimeVector
!</inputoutput>

!</subroutine>
    ! local variables
    real(DP), dimension(:), pointer :: p_Dsource

    ! Make sure we can store the timestep data.
    if ((isubvector .lt. 1) .or. (isubvector .gt. rspaceTimeVector%NEQtime)) then
      call output_line('Invalid timestep number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_setTimestepData')
      call sys_halt()
    end if
    
    if (rvector%NEQ .ne. rspaceTimeVector%NEQ) then
      call output_line('Vector size invalid!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_setTimestepData')
      call sys_halt()
    end if

    ! Save the vector data. If necessary, new memory is allocated -- as the
    ! default value of the handle in the handle list is ST_NOHANDLE.
    !
    ! Don't use storage_copy, as this might give errors in case the array
    ! behind the handle is longer than the vector!
    call lsysbl_getbase_double (rvector,p_Dsource)
    !CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(isubvector),p_Ddest)
    !CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
    call exstor_setdata_double (rspaceTimeVector%p_IdataHandleList(isubvector),p_Dsource)

    ! After a setTimestepData, the scale factor is 1.0.
    rspaceTimeVector%p_Dscale(isubvector) = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_getTimestepData (rspaceTimeVector, isubvector, rvector)

!<description>
  ! Restores the data of timestep itimestep into the vector rvector.
!</desctiprion>

!<input>
  ! Number of the subvector that corresponds to rvector. >= 1, <= NEQtime from
  ! the initialisation of rspaceTimeVector.
  integer, intent(IN) :: isubvector
  
  ! Space-time vector structure where to save the data.
  type(t_spacetimeVector), intent(IN) :: rspaceTimeVector
!</input>

!<inputoutput>
  ! Vector with data that should receive the data.
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddest

    ! Make sure we can store the timestep data.
    if ((isubvector .lt. 1) .or. (isubvector .gt. rspaceTimeVector%NEQtime)) then
      call output_line('Invalid timestep number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      call sys_halt()
    end if
    
    if (rvector%NEQ .ne. rspaceTimeVector%NEQ) then
      call output_line('Vector size invalid!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      call sys_halt()
    end if
    
    if (rspaceTimeVector%p_Dscale(isubvector) .eq. 0.0_DP) then
     ! The vector is a zero vector
      call lsysbl_clearVector (rvector)
      return
    end if

    ! Get the vector data. 
    !
    ! Don't use storage_copy, as this might give errors in case the array
    ! behind the handle is longer than the vector!
    !CALL storage_getbase_double (rspaceTimeVector%p_IdataHandleList(isubvector),&
    !    p_Dsource)
    call lsysbl_getbase_double (rvector,p_Ddest)
    !CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
    call exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(isubvector),p_Ddest)
    
    ! Scale the vector?
    if (rspaceTimeVector%p_Dscale(isubvector) .ne. 1.0_DP) then
      call lalg_scaleVectorDble (p_Ddest,rspaceTimeVector%p_Dscale(isubvector))
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_getTimestepDataByTime (rspaceTimeVector, dtimestamp, rvector)

!<description>
  ! Restores the data of a timestep into the vector rvector. dtimestamp is a time
  ! stamp in the range 0.0 .. 1.0, where 0.0 corresponds to the 1st subvector
  ! and 1.0 to the last subvector in the space time vector. If dtimestamp
  ! specifies a time stamp between two stored vectors, quadratic interpolation
  ! is used to calculate rvector.
!</desctiprion>

!<input>
  ! Time stamp of the vector whose data should be retrieved.
  real(DP), intent(IN) :: dtimestamp
  
  ! Space-time vector structure where to save the data.
  type(t_spacetimeVector), intent(IN) :: rspaceTimeVector
!</input>

!<inputoutput>
  ! Vector with data that should receive the data.
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddest
    integer :: itimestep1,itimestep2,itimestep3
    real(DP) :: dreltime,dabstime
    integer :: i
    real(DP) :: dscal1,dscal2,dscal3
    real(DP), dimension(:,:), allocatable :: p_Dsource

    ! Make sure we can store the timestep data.
    if ((dtimestamp .lt. 0.0_DP) .or. (dtimestamp .gt. 1.0_DP)) then
      call output_line('Invalid time stamp!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepDataByTime')
      call sys_halt()
    end if
    
    if (rvector%NEQ .ne. rspaceTimeVector%NEQ) then
      call output_line('Vector size invalid!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_getTimestepData')
      call sys_halt()
    end if
    
    ! Get the time step which is closest to the time stamp.
    ! Rescale dtimestamp to the interval [1..NEQtime].
    dabstime = dtimestamp*real(rspaceTimeVector%NEQtime-1,DP)+1.0_DP
    itimestep2 = int(dabstime + 0.5_DP)
    
    if (dabstime .eq. real(itimestep2,DP)) then
      ! Nice coincidence, we have exactly timestep itimestep2. Ok, then we 
      ! can call the routine to get that timestep; this saves us some
      ! time as the interpolation can be omitted.
      call sptivec_getTimestepData (rspaceTimeVector, itimestep2, rvector)
      return
    end if
    
    if (rspaceTimeVector%NEQtime .eq. 2) then
      ! Special case: only one timestep!
      itimestep1 = 0
      itimestep2 = 0
      itimestep3 = 1
    else
      ! Is this the first or the last timestep?
      if (itimestep2 .eq. 1) then
        ! First timestep. Interpolate between timesteps 0,1 and 2, evaluate 
        ! near timestep 0.
        itimestep1 = 0
        itimestep2 = 1
        itimestep3 = 2
      else if (itimestep2 .eq. rspaceTimeVector%NEQtime) then
        ! Last timestep. Interpolate between timesteps n-2,n-1 and n, evaluate 
        ! near timestep n.
        itimestep1 = rspaceTimeVector%NEQtime-2
        itimestep2 = rspaceTimeVector%NEQtime-1
        itimestep3 = rspaceTimeVector%NEQtime
      else
        ! Somewhere in the inner. Get the number of the previous and next timestep
        itimestep1 = itimestep2-1
        itimestep3 = itimestep2+1
      end if
    end if

    ! Calculate the 'relative' time in the interval [-1,1], where -1 corresponds
    ! to timestep itimestep1, 0 to itimestep2 and +1 to itimestep3.
    ! This will be used to evaluate the quadratic polynomial.
    dreltime = dabstime-1.0_DP-real(itimestep2,DP)
    
    ! Get the vector data of the three timesteps
    allocate(p_Dsource(rspaceTimeVector%NEQ,3))
    
    call exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(itimestep1),&
        p_Dsource(:,1))
    call exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(itimestep2),&
        p_Dsource(:,2))
    call exstor_getdata_double (rspaceTimeVector%p_IdataHandleList(itimestep3),&
        p_Dsource(:,3))

    call lsysbl_getbase_double (rvector,p_Ddest)

    ! Calculate the quadratic interpolation of the three arrays.
    dscal1 = rspaceTimeVector%p_Dscale(itimestep1)
    dscal2 = rspaceTimeVector%p_Dscale(itimestep2)
    dscal3 = rspaceTimeVector%p_Dscale(itimestep3)
    if (rspaceTimeVector%NEQtime .eq. 2) then
      ! Special case: only 1 timestep. Linear interpolation. dreltime is in [0..1]!
      do i=1,size(p_Ddest)
        p_Ddest(i) = (1.0_DP-dreltime) * dscal2*p_Dsource(i,2) + &
                     dreltime*dscal3*p_Dsource(i,3)
      end do
    else
      ! Quadratic interpolation
      do i=1,size(p_Ddest)
        call mprim_quadraticInterpolation (dreltime,&
            dscal1*p_Dsource(i,1),dscal2*p_Dsource(i,2),dscal3*p_Dsource(i,3),p_Ddest(i))
      end do
    end if
    
    deallocate(p_Dsource)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_vectorLinearComb (rx,ry,cx,cy)

!<description>
  ! Performs a linear combination of space-time vectors: ry = cx * rx  +  cy * ry
!</desctiprion>

!<input>
  ! First source vector
  type(t_spacetimeVector), intent(IN)   :: rx
  
  ! Scaling factor for Dx
  real(DP), intent(IN)               :: cx

  ! Scaling factor for Dy
  real(DP), intent(IN)               :: cy
!</input>

!<inputoutput>
  ! Second source vector; also receives the result
  type(t_spacetimeVector), intent(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

    integer :: i
    integer(PREC_VECIDX), dimension(1) :: Isize
    type(t_vectorBlock) :: rxBlock,ryBlock

    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    
    if (rx%NEQ .ne. ry%NEQ) then
      call output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_vectorLinearComb')
      call sys_halt()
    end if

    if (rx%NEQtime .ne. ry%NEQtime) then
      call output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_vectorLinearComb')
      call sys_halt()
    end if
    
    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    call lsysbl_createVecBlockDirect (rxBlock,Isize,.false.)
    call lsysbl_createVecBlockDirect (ryBlock,Isize,.false.)

    ! DEBUG!!!
    call lsysbl_getbase_double (rxBlock,p_Dx)
    call lsysbl_getbase_double (ryBlock,p_Dy)

    ! Loop through the substeps, load the data in, perform the linear combination
    ! and write out again.
    do i=1,rx%NEQtime
      call sptivec_getTimestepData (rx, i, rxBlock)
      call sptivec_getTimestepData (ry, i, ryBlock)
      
      call lsysbl_vectorLinearComb (rxBlock,ryBlock,cx,cy)

      call sptivec_setTimestepData (ry, i, ryBlock)
    end do

    ! Release temp memory    
    call lsysbl_releaseVector (ryBlock)
    call lsysbl_releaseVector (rxBlock)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_copyVector (rx,ry)

!<description>
  ! Copys a vector: ry := rx.
!</desctiprion>

!<input>
  ! Source vector
  type(t_spacetimeVector), intent(IN)   :: rx
!</input>

!<inputoutput>
  ! Destination vector
  type(t_spacetimeVector), intent(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

    integer :: i
    integer(PREC_VECIDX), dimension(1) :: Isize
    type(t_vectorBlock) :: rxBlock,ryBlock
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    
    if ((ry%NEQ .eq. 0) .and. (ry%NEQtime .eq. 0)) then
      ! Destination vector does not exist. Create it.
      call sptivec_initVector (ry,rx%NEQtime,rx%p_rblockDiscretisation)
    end if
    
    if (rx%NEQ .ne. ry%NEQ) then
      call output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_copyVector')
      call sys_halt()
    end if

    if (rx%NEQtime .ne. ry%NEQtime) then
      call output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_copyVector')
      call sys_halt()
    end if

    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    call lsysbl_createVecBlockDirect (rxBlock,Isize,.false.)
    call lsysbl_createVecBlockDirect (ryBlock,Isize,.false.)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rxBlock,p_Dx)
    call lsysbl_getbase_double (ryBlock,p_Dy)

    ! Loop through the substeps, load the data in, perform the linear combination
    ! and write out again.
    do i=1,rx%NEQtime
      
      if (rx%p_Dscale(i) .eq. 0.0_DP) then
        ry%p_Dscale(i) = 0.0_DP
      else
        call sptivec_getTimestepData (rx, i, rxBlock)

        call lsysbl_copyVector (rxBlock,ryBlock)

        call sptivec_setTimestepData (ry, i, ryBlock)
      end if

    end do

    ! Release temp memory    
    call lsysbl_releaseVector (ryBlock)
    call lsysbl_releaseVector (rxBlock)
      
  end subroutine

  ! ***************************************************************************

!<function>
  
  real(DP) function sptivec_scalarProduct (rx, ry)
  
!<description>
  ! Calculates a scalar product of two block vectors.
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>

!<input>
  ! First source vector
  type(t_spacetimeVector), intent(IN)   :: rx

  ! Second source vector
  type(t_spacetimeVector), intent(IN)   :: ry
!</input>

!<result>
  ! The scalar product (rx,ry) of the two block vectors.
!</result>

!</function>

    integer :: i
    integer(PREC_VECIDX), dimension(1) :: Isize
    real(DP) :: dres
    type(t_vectorBlock) :: rxBlock,ryBlock
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    
    if (rx%NEQ .ne. ry%NEQ) then
      call output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      call sys_halt()
    end if

    if (rx%NEQtime .ne. ry%NEQtime) then
      call output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      call sys_halt()
    end if

    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    call lsysbl_createVecBlockDirect (rxBlock,Isize,.false.)
    call lsysbl_createVecBlockDirect (ryBlock,Isize,.false.)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rxBlock,p_Dx)
    call lsysbl_getbase_double (ryBlock,p_Dy)

    ! Loop through the substeps, load the data in, perform the scalar product.
    dres = 0.0_DP
    do i=1,rx%NEQtime
      
      if ((rx%p_Dscale(i) .ne. 0.0_DP) .and. (ry%p_Dscale(i) .ne. 0.0_DP)) then
        call sptivec_getTimestepData (rx, i, rxBlock)
        call sptivec_getTimestepData (ry, i, ryBlock)
        
        dres = dres + lsysbl_scalarProduct (rxBlock,ryBlock)
      end if

    end do

    ! Release temp memory    
    call lsysbl_releaseVector (ryBlock)
    call lsysbl_releaseVector (rxBlock)
      
    sptivec_scalarProduct = dres
      
  end function

  ! ***************************************************************************

!<function>
  
  real(DP) function sptivec_scalarProductWeighted (rx, ry, Dweights)
  
!<description>
  ! Calculates a weighted scalar product of two block vectors.
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>

!<input>
  ! First source vector
  type(t_spacetimeVector), intent(IN)   :: rx

  ! Second source vector
  type(t_spacetimeVector), intent(IN)   :: ry
  
  ! An array with weights for each component of a solution
  ! vector in a timestep. Dweights(i) is multiplied to
  ! the i'th solution component of each timestep.
  real(DP), dimension(:), intent(IN) :: Dweights
!</input>

!<result>
  ! The scalar product (rx,ry) of the two block vectors.
!</result>

!</function>

    integer :: i,irow
    integer(PREC_VECIDX), dimension(1) :: Isize
    real(DP) :: dres,a
    type(t_vectorBlock) :: rxBlock,ryBlock
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    
    if (rx%NEQ .ne. ry%NEQ) then
      call output_line('Space-time vectors have different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      call sys_halt()
    end if

    if (rx%NEQtime .ne. ry%NEQtime) then
      call output_line('Space-time vectors have different number of timesteps!',&
          OU_CLASS_ERROR,OU_MODE_STD,'sptivec_scalarProduct')
      call sys_halt()
    end if

    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    call lsysbl_createVecBlockDirect (rxBlock,Isize,.false.)
    call lsysbl_createVecBlockDirect (ryBlock,Isize,.false.)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rxBlock,p_Dx)
    call lsysbl_getbase_double (ryBlock,p_Dy)

    ! Loop through the substeps, load the data in, perform the scalar product.
    dres = 0.0_DP
    do i=1,rx%NEQtime
      
      if ((rx%p_Dscale(i) .ne. 0.0_DP) .and. (ry%p_Dscale(i) .ne. 0.0_DP)) then
        call sptivec_getTimestepData (rx, i, rxBlock)
        call sptivec_getTimestepData (ry, i, ryBlock)
        
        ! Calculate a weighted scalar product using the equation weights.
        a = 0.0_DP
        do irow = 1,rxBlock%nblocks
          a = a + Dweights(irow) * &
              lsyssc_scalarProduct(rxBlock%RvectorBlock(irow),&
                                   ryBlock%RvectorBlock(irow))
        end do
        
        dres = dres + a
      end if

    end do

    ! Release temp memory    
    call lsysbl_releaseVector (ryBlock)
    call lsysbl_releaseVector (rxBlock)
      
    sptivec_scalarProductWeighted = dres
      
  end function

  ! ***************************************************************************

!<function>

  real(DP) function sptivec_vectorNorm (rx,cnorm)

!<description>
  ! Calculates the norm of the vector rx.
!</desctiprion>

!<input>
  ! Source vector
  type(t_spacetimeVector), intent(IN)   :: rx

  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm
!</input>

!<result>
  ! Norm of the vector.
!</result>
  
!</function>

    integer(PREC_VECIDX), dimension(1) :: Isize
    type(t_vectorBlock) :: rxBlock
    real(DP) :: dnorm
    integer :: i

    Isize(1) = rx%NEQ

    ! Allocate a 'little bit' of memory for the subvectors
    call lsysbl_createVecBlockDirect (rxBlock,Isize,.false.)
    
    dnorm = 0.0_DP
    
    ! Loop through the substeps, load the data in, sum up to the norm.
    do i=1,rx%NEQtime
      if (rx%p_Dscale(i) .ne. 0.0_DP) then
        call sptivec_getTimestepData (rx, i, rxBlock)
        
        select case (cnorm)
        case (LINALG_NORML2)
          dnorm = dnorm + lsysbl_vectorNorm (rxBlock,cnorm)**2
        case DEFAULT
          dnorm = dnorm + lsysbl_vectorNorm (rxBlock,cnorm)
        end select
      end if
    end do

    ! Release temp memory    
    call lsysbl_releaseVector (rxBlock)
    
    ! Calculate the actual norm.
    select case (cnorm)
    case (LINALG_NORML1)
      sptivec_vectorNorm = dnorm / (rx%NEQtime)
    case (LINALG_NORML2)
      sptivec_vectorNorm = sqrt(dnorm / (rx%NEQtime))
    case DEFAULT
      sptivec_vectorNorm = dnorm
    end select

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_clearVector (rx)

!<description>
  ! Initialises a vector with zero.
!</desctiprion>

!<inputoutput>
  ! Vector to be cleared.
  type(t_spacetimeVector), intent(INOUT)   :: rx
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
    call lalg_clearVectorDble (rx%p_Dscale(:))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_scaleVector (rx,dscale)

!<description>
  ! Scales a vector by dscale.
!</desctiprion>

!<input>
  ! Scaling factor.
  real(DP), intent(IN) :: dscale
!</input>

!<inputoutput>
  ! Vector to be scaled.
  type(t_spacetimeVector), intent(INOUT)   :: rx
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
    call lalg_scaleVectorDble(rx%p_Dscale(:),dscale)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_convertSupervecToVector (rxsuper, rx)

!<description>
  ! Converts a space-time coupled vector into a usual block vector.
  ! The blocks in rx correspond to the timesteps in rxsuper.
  !
  ! WARNING: This may be memory intensive!
!</description>

!<input>
  ! Space time vector to be converted.
  type(t_spacetimeVector), intent(IN)   :: rxsuper
!</input>

!<inputoutput>
  ! Destination block vector that should receive the result.
  type(t_vectorBlock), intent(INOUT) :: rx
!</inputoutput>

!</subroutine>

    integer(PREC_VECIDX), dimension(:), allocatable :: Isize
    type(t_vectorBlock) :: rvectorTmp
    integer :: i
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    
    ! Create a complete new vector
    if (rx%NEQ .ne. 0) call lsysbl_releaseVector (rx)
  
    ! Create a vector in the correct size
    allocate(Isize(rxsuper%NEQtime))
    Isize(:) = rxsuper%NEQ
    call lsysbl_createVecBlockDirect (rx,Isize,.false.)

    ! Create a 1-block temp vector for the data    
    call lsysbl_createVecBlockDirect (rvectorTmp,Isize(1:1),.false.)
    call lsysbl_getbase_double (rvectorTmp,p_Ddata1)
    
    ! Load the subvectors and write them to the global vector.
    do i=1,rxsuper%NEQtime
      call sptivec_getTimestepData (rxsuper, i, rvectorTmp)
      
      call lsyssc_getbase_double (rx%RvectorBlock(i),p_Ddata2)
      call lalg_copyVectorDble (p_Ddata1,p_Ddata2)
    end do
    
    ! Release the temp vector
    call lsysbl_releaseVector (rvectorTmp)
    
    deallocate(Isize)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_convertVectorToSupervec (rx,rxsuper)

!<description>
  ! Converts a usual block vector into a space-time coupled vector.
  ! The blocks in rx correspond to the timesteps in rxsuper.
!</description>

!<input>
  ! Block vector to be converted.
  type(t_vectorBlock), intent(IN)   :: rx
!</input>

!<inputoutput>
  ! Destination space-time vector that should receive the result.
  type(t_spacetimeVector), intent(INOUT) :: rxsuper
!</inputoutput>

!</subroutine>

    integer(PREC_VECIDX), dimension(1) :: Isize
    type(t_vectorBlock) :: rvectorTmp
    integer :: i
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    
    ! Create a 1-block temp vector for the data    
    Isize(1) = rxsuper%NEQ
    call lsysbl_createVecBlockDirect (rvectorTmp,Isize(1:1),.false.)
    call lsysbl_getbase_double (rvectorTmp,p_Ddata2)
    
    ! Load the subvectors and write them to the global vector.
    do i=1,rxsuper%NEQtime
      call lsyssc_getbase_double (rx%RvectorBlock(i),p_Ddata1)
      call lalg_copyVectorDble (p_Ddata1,p_Ddata2)

      call sptivec_setTimestepData (rxsuper, i, rvectorTmp)
    end do
    
    ! Release the temp vector
    call lsysbl_releaseVector (rvectorTmp)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_saveToCollection (rx,rcollection,sname,ilevel,ssection)

!<description>
  ! Saves a space-time vector to a collection.
!</description>

!<input>
  ! Space-time vector to be saved.
  type(t_spacetimeVector), intent(IN)   :: rx
  
  ! Name that the vector should be given in the collection.
  character(LEN=*), intent(IN) :: sname

  ! OPTIONAL: Name of the section in the collection where to save the vector to.
  character(LEN=*), intent(IN), optional :: ssection
  
  ! OPTIONAL: Level in the collection where to save the vector to.
  integer, intent(IN), optional :: ilevel
!</input>

!<inputoutput>
  ! Collection structure. The vector is saved to this.
  type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!</subroutine>

    ! Save the content of the structure    
    call collct_setvalue_int (rcollection,trim(sname)//'_NEQ',&
        rx%NEQ,.true.,ilevel,ssection)
  
    call collct_setvalue_int (rcollection,trim(sname)//'_NTST',&
        rx%NEQtime,.true.,ilevel,ssection)

    if (rx%NEQtime .ne. 0) then
      call collct_setvalue_intarr (rcollection,trim(sname)//'_NTST',&
          rx%p_IdataHandleList,.true.,ilevel,ssection)
      
      call collct_setvalue_realarr (rcollection,trim(sname)//'_SCALE',&
          rx%p_Dscale,.true.,ilevel,ssection)

      ! Otherwise the pointers are NULL()!
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_restoreFromCollection (rx,rcollection,sname,ilevel,ssection)

!<description>
  ! Restores a space-time vector from a collection.
  !
  ! Note that this creates a copy of a previous space-time vector. This vector
  ! must be released with releaseVector to prevent memory leaks
!</description>

!<input>
  ! Collection structure where to restore rx from.
  type(t_collection), intent(INOUT) :: rcollection

  ! Name that the vector should be given in the collection.
  character(LEN=*), intent(IN) :: sname

  ! OPTIONAL: Name of the section in the collection where to save the vector to.
  character(LEN=*), intent(IN), optional :: ssection
  
  ! OPTIONAL: Level in the collection where to save the vector to.
  integer, intent(IN), optional :: ilevel
!</input>

!<inputoutput>
  ! Space-time vector that receives the data.
  type(t_spacetimeVector), intent(OUT)   :: rx
!</inputoutput>

!</subroutine>
    
    ! The vector is a copy of another one, so the releaseVector routine
    ! will not release the content.
    rx%bisCopy = .true.

    ! Get the content of the structure    
    rx%NEQ = collct_getvalue_int (rcollection,trim(sname)//'_NEQ',&
        ilevel,ssection)
  
    rx%NEQtime = collct_getvalue_int (rcollection,trim(sname)//'_NTST',&
        ilevel,ssection)
        
    if (rx%NEQtime .ne. 0) then
      ! For the handle list, we need to allocate some memory...
      allocate(rx%p_IdataHandleList(rx%NEQtime))

      call collct_getvalue_intarr (rcollection,trim(sname)//'_NTST',&
          rx%p_IdataHandleList,ilevel,ssection)

      call collct_getvalue_realarr (rcollection,trim(sname)//'_SCALE',&
          rx%p_Dscale,ilevel,ssection)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_removeFromCollection (rcollection,sname,ilevel,ssection)

!<description>
  ! Removes a space-time vector from a collection.
!</description>

!<input>
  ! Name that the vector should be given in the collection.
  character(LEN=*), intent(IN) :: sname

  ! OPTIONAL: Name of the section in the collection where to save the vector to.
  character(LEN=*), intent(IN), optional :: ssection
  
  ! OPTIONAL: Level in the collection where to save the vector to.
  integer, intent(IN), optional :: ilevel
!</input>

!<inputoutput>
  ! Collection structure where to remove sname from.
  type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!</subroutine>

    ! Remove all entries
    call collct_deleteValue (rcollection,trim(sname)//'_DIR',&
        ilevel,ssection)

    call collct_deleteValue (rcollection,trim(sname)//'_SCALE',&
        ilevel,ssection)
  
    call collct_deleteValue (rcollection,trim(sname)//'_NEQ',&
        ilevel,ssection)
  
    call collct_deleteValue (rcollection,trim(sname)//'_NTST',&
        ilevel,ssection)
        
    call collct_deleteValue (rcollection,trim(sname)//'_NTST',&
        ilevel,ssection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_loadFromFileSequence (rx,sfilename,istart,iend,idelta,&
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
  character(LEN=*), intent(IN) :: sfilename

  ! Number of the first file to be read in
  integer, intent(IN) :: istart
  
  ! Number of the last file to be read in
  integer, intent(IN) :: iend
  
  ! Delta parameter that specifies how to increase the filename suffix.
  ! Standard is =1.
  integer, intent(IN) :: idelta
  
  ! Whether to read formatted or unformatted data from disc.
  logical, intent(IN) :: bformatted

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
  logical, optional :: brepeatLast
  
  ! OPTIONAL: A block discretisation structure that defines the shape of the spatial
  ! vectors. If not specified, the vectors will be read in as pure data vectors
  ! without a discretisation attached.
  type(t_blockDiscretisation), intent(IN), target, optional :: rblockDiscretisation

!</input>

!<inputoutput>
  ! Space-time vector where store data to.
  ! If empty, the vector is created from the scratch.
  ! If the vector is not empty, it must be large enough to hold all data.
  ! The content is ovrwritten.
  type(t_spaceTimeVector), intent(INOUT) :: rx
!</inputoutput>

!</subroutine>

    ! Local variables
    type(t_vectorBlock) :: rvector
    type(t_vectorScalar) :: rvectorScalar
    character(SYS_STRLEN) :: sfile,sarray
    integer :: i,ilast
    logical :: bexists,brepeat
    
    brepeat = .false.
    if (present(brepeatLast)) brepeat = brepeatLast
    
    ilast = 0

    ! Loop over the files
    do i=istart,iend
    
      ! Form the filename
      write(sfile,sfilename) i*idelta
      
      ! Is the file there?
      inquire(file=trim(sfile), exist=bexists)
      
      if (bexists) then
        ! Remember this solution as the last available one
        ilast = i
      
        if (.not. present(rblockDiscretisation)) then
          ! Read the file into rvector. The first read command creates rvector
          ! in the correct size.
          call vecio_readVectorHR (rvectorScalar, sarray, .false.,&
            0, sfile, bformatted)

          if ((i .eq. istart) .and. (rx%NEQtime .eq. 0)) then
            ! At the first file, create a space-time vector holding the data.
            call sptivec_initVectorPlain (rx,rvectorScalar%NEQ,iend-istart+1)
          end if         

          ! Save the data
          call exstor_setdata_storage (rx%p_IdataHandleList(1+i),rvectorScalar%h_Ddata)
    
        else
  
          call vecio_readBlockVectorHR (rvector, sarray, .false.,&
            0, sfile, bformatted)

          if (i .eq. istart) then
            ! At the first file, create a space-time vector holding the data.
            call sptivec_initVector (rx,1+iend-istart,rblockDiscretisation)
          end if         
          ! Save the data
          call exstor_setdata_storage (rx%p_IdataHandleList(1+i),rvector%h_Ddata)
          
        end if
        
      else
        if (i .eq. istart) then
          ! The first file must exist!
          call output_line ('The first file must exist!', &
              ssubroutine='sptivec_loadFromFileSequence')
        end if         

        if (brepeat) then
          call output_line ('Warning: Unable to load file "'//trim(sfile) &
              //'". Repeating last solution!', &
              ssubroutine='sptivec_loadFromFileSequence')
        
          ! Copy the data from the last known solution to the current one.
          call exstor_copy (&
              rx%p_IdataHandleList(1+ilast),&
              rx%p_IdataHandleList(1+i))
        else
          call output_line ('Warning: Unable to load file "'//trim(sfile) &
              //'". Assuming zero!', ssubroutine='sptivec_loadFromFileSequence')
        
          ! Clear that array. Zero solution.
          call exstor_clear (rx%p_IdataHandleList(1+i))
        end if
      end if
    
    end do
    
    ! The vector is scaled by 1.0.
    ! rx%p_Dscale(:) = 1.0_DP
    call lalg_setVectorDble(rx%p_Dscale(:),1.0_DP)
    
    ! Remove the temp vector
    if (present(rblockDiscretisation)) then
      call lsysbl_releaseVector (rvector)
    else
      call lsyssc_releaseVector (rvectorScalar)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_saveToFileSequence (rx,sfilename,bformatted,&
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
  type(t_spaceTimeVector), intent(IN) :: rx

  ! Filename pattern + path where to form a filename from.
  character(LEN=*), intent(IN) :: sfilename

  ! Whether to read formatted or unformatted data from disc.
  logical, intent(IN) :: bformatted
  
  ! OPTIONAL: Format string that is used for exporting data to files.
  ! E.g. '(E20.10)'.
  character(LEN=SYS_STRLEN), intent(IN), optional :: sformat
!</input>

!<inputoutput>
  ! Temporary vector. This vector must prescribe the block structure of the
  ! subvectors in the space-time vector. If not specified, the data is written
  ! out without a block structure.
  type(t_vectorBlock), intent(INOUT), target, optional :: rtempVector
!</inputoutput>

!</subroutine>

    ! Local variables
    character(SYS_STRLEN) :: sfile
    integer :: i
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx
    type(t_vectorBlock), pointer :: p_rx
    integer(PREC_VECIDX), dimension(1) :: Isize
    
    if (present(rtempVector)) then
      p_rx => rtempVector
    else
      ! Create a 1-block temp vector for the data    
      allocate(p_rx)
      Isize(1) = rx%NEQ
      call lsysbl_createVecBlockDirect (p_rx,Isize(1:1),.false.)
    end if

    ! DEBUG!!!
    call lsysbl_getbase_double (p_rx,p_Dx)

    ! Loop over the files
    do i=0,rx%NEQtime-1
    
      ! Get the data from the space-time vector.
      call sptivec_getTimestepData (rx, 1+i, p_rx)
      
      ! Form the filename
      write(sfile,sfilename) i
      
      ! Save that to disc.
      if (.not. bformatted) then
        call vecio_writeBlockVectorHR (p_rx, 'vector', .true.,&
                                      0, sfile)
      else if (.not. present(sformat)) then
          call vecio_writeBlockVectorHR (p_rx, 'vector', .true.,&
                                        0, sfile, '(E24.16)')
      else
        call vecio_writeBlockVectorHR (p_rx, 'vector', .true.,&
                                       0, sfile, sformat)
      end if
    
    end do
    
    if (.not. present(rtempVector)) then
      call lsysbl_releaseVector (p_rx)
      deallocate(p_rx)
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptivec_setConstant (rx,dvalue)

!<description>
  ! Sets the whole space-time vector to a constant value.
!</description>

!<input>
  ! Value which should be written to the whole space-time vector.
  real(DP), intent(IN) :: dvalue
!</input>

!<inputoutput>
  ! Space-time vector to modify
  type(t_spaceTimeVector), intent(INOUT) :: rx
!</inputoutput>

!</subroutine>

    ! Local variables
    real(DP), dimension(:), allocatable :: p_Ddata
    integer :: i
    
    ! Allocate memory for intermediate values
    allocate(p_Ddata(rx%NEQ))

    ! Loop over the files
    do i=1,rx%NEQtime
    
      ! Get the data and set to a defined value.
      call exstor_getdata_double (rx%p_IdataHandleList(i),p_Ddata)
      p_Ddata(:) = dvalue
      call exstor_setdata_double (rx%p_IdataHandleList(i),p_Ddata)

      rx%p_Dscale(i) = 1.0_DP

    end do
    
    deallocate(p_Ddata)
    
  end subroutine

end module
