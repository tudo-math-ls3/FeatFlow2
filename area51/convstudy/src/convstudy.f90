!##############################################################################
!# ****************************************************************************
!# <name> convstudy </name>
!# ****************************************************************************
!#
!# <purpose>
!# This programs can be used to perform convergence studies on a
!# sequence of successively refined grids. This program accepts two
!# different GMV file representing finite element solutions on a
!# coarse and a fine grid and computes the different error norms for
!# each common variable present in both data sets.
!#
!# </purpose>
!##############################################################################

program convstudy

  use basicgeometry
  use collection
  use convstudy_callback
  use derivatives
  use domainintegration
  use element
  use feevaluation
  use fsystem
  use genoutput
  use linearsystemscalar
  use linearsystemblock
  use pprocerror
  use spatialdiscretisation
  use storage
  use triangulation
  use ucd

  implicit none

!<types>

!<typeblock>
  
  ! Data set structure. This structure holds all data read in from the GMV file

  type t_dataSet

    ! Triangulation structure
    type(t_triangulation) :: rtriangulation

    ! Discretisation structure
    type(t_spatialDiscretisation) :: rdiscretisation

    ! Scalar solution vectors
    type(t_vectorScalar), dimension(:), pointer :: Rvector

    ! Names of scalar solution vectors
    character(LEN=SYS_STRLEN), dimension(:), pointer :: SvarNames

  end type t_dataSet

!</typeblock>

!</types>

  ! local variables
  type(t_dataSet) :: rdataSet, rdataSetRef

  !*****************************************************************************

  ! Initialize Feat2 subsystem
  call system_init()

  ! Initialize storage subsystem
  call storage_init(500, 100)

  ! Read in data sets
  call convst_readDataSets(rdataSet, rdataSetRef)

  ! Compute error norms
  call convst_computeErrorNorms(rdataSet, rdataSetRef)

  ! Release data sets
  call convst_releaseDataSet(rdataSet)
  call convst_releaseDataSet(rdataSetRef)

  ! Release storage
  call storage_done()
  call output_lbrk()

contains

  !*****************************************************************************

!<subroutine>

  subroutine convst_readDatasets(rdataSet, rdataSetRef)

!<description>
    ! This subroutine reads in the data sets from the GMV file given
    ! as command line argument. If no command line arguments are given
    ! an error is thrown.
!</description>

!<output>
    ! Data sets to be filled with data from GMV files
    type(t_dataSet), intent(out) :: rdataSet, rdataSetRef
!</output>
!</subroutine>
  
    ! local variable
    character(LEN=SYS_STRLEN) :: cbuffer
    integer :: iarg, narg

    iarg = 1; narg = command_argument_count()

    cmdarg: do
      ! Retrieve next command line argument
      call get_command_argument(iarg,cbuffer)

      if (trim(adjustl(cbuffer)) .eq. '--ucdfile') then
        
        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call convst_readDataset(cbuffer, rdataSet)

      elseif (trim(adjustl(cbuffer)) .eq. '--ucdreffile') then

        iarg = iarg+1
        call get_command_argument(iarg,cbuffer)
        call convst_readDataset(cbuffer, rdataSetRef)

      else
        iarg = iarg+1
        if (iarg .ge. narg) exit cmdarg
      end if
    end do cmdarg
  
  end subroutine convst_readDatasets

  !*****************************************************************************

!<subroutine>

    subroutine convst_readDataset(sfilename, rdataSet)

!<description>
    ! This subroutine reads in a single data sets from the GMV file.
!</description>

!<input>
  ! Name of the GMV file.
  character(LEN=*), intent(in) :: sfilename
!</input>

!<output>
    ! Data set to be filled with data from GMV files
    type(t_dataSet), intent(out) :: rdataSet
!</output>
!</subroutine>
    
    ! local variable
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: i,nlength

    ! Read GMV file into UCD export structure
    call ucd_readGMV(sfilename, rexport, rdataSet%rtriangulation)

    ! What spatial dimension are we?
    select case(rdataSet%rtriangulation%ndim)

    case (NDIM1D)
      ! Create spatial discretisation in 1D
      call spdiscr_initDiscr_simple(rdataSet%rdiscretisation, EL_E001_1D,&
          SPDISC_CUB_AUTOMATIC, rdataSet%rtriangulation)

    case (NDIM2D)
      ! Create spatial discretisation in 2D
      call spdiscr_initDiscr_triquad(rdataSet%rdiscretisation, EL_E001, EL_E011,&
          SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, rdataSet%rtriangulation)

    case (NDIM3D)
!!$      ! Create spatial discretisation in 3D
!!$      call spdiscr_initDiscr_simple(rdataSet%rdiscretisation, EL_E001_3D,&
!!$          SPDISC_CUB_AUTOMATIC, rdataSet%rtriangulation)
      
      ! Create spatial discretisation in 3D
      call spdiscr_initDiscr_simple(rdataSet%rdiscretisation, EL_E011_3D,&
          SPDISC_CUB_AUTOMATIC, rdataSet%rtriangulation)
    end select

    ! Allocate memory for scalar solution vectors
    allocate(rdataSet%Rvector(rexport%nvariables))
    allocate(rdataSet%SvarNames(rexport%nvariables))
    
    ! Copy scalar solution vectors
    do i = 1, rexport%nvariables
      
      call lsyssc_createVector(rdataSet%rdiscretisation, rdataSet%Rvector(i))
      call lsyssc_getbase_double(rdataSet%Rvector(i), p_Ddata)

      rdataSet%SvarNames(i) = rexport%p_SvariableNames(i)
      call ucd_getVariable(rexport, rdataSet%SvarNames(i), p_Ddata, nlength)

      if (nlength .ne. size(p_Ddata)) then
        call output_line('Invalid vector size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convst_readDataset')
        call sys_halt()
      end if

    end do

    ! Release UCD export structure
    call ucd_release(rexport)

  end subroutine convst_readDataset

  !*****************************************************************************

!<subroutine>

  subroutine convst_computeErrorNorms(rdataSet, rdataSetRef)

!<description>
    ! This subroutine computes the norms of the L1-, L2-, and H1-errors
    ! of each common component of the given data sets
!</description>

!<input>
    ! Data sets for which to compute the error norms
    type(t_dataSet), intent(in), target :: rdataSet, rdataSetRef
!</input>

    ! local variables
    type(t_collection) :: rcollection
    type(t_vectorBlock), target :: rvector
    real(DP) :: derrorL1, derrorL2, derrorH1
    integer :: i,j 

    do i = 1, size(rdataSetRef%Rvector)
      do j = 1, size(rdataSet%Rvector)

        if (trim(adjustl(sys_upcase(rdataSetRef%SvarNames(i)))) .eq.&
            trim(adjustl(sys_upcase(rdataSet%SvarNames(j))))) then

          ! Create temporal 1-block vectors
          call lsysbl_createVecFromScalar(rdataSetRef%Rvector(i), rvector)

          ! Prepare collection
          rcollection%IquickAccess(1) = PPERR_L1ERROR
          rcollection%p_rvectorQuickAccess1 => rvector
          
          ! Compute L1-error
          call pperr_scalar(PPERR_L1ERROR, derrorL1,&
              rdataSet%Rvector(j), convst_refFunction, rcollection)
          
          ! Compute L2-error
          call pperr_scalar(PPERR_L2ERROR, derrorL2,&
              rdataSet%Rvector(j), convst_refFunction, rcollection)

          ! Compute H1-error
          call pperr_scalar(PPERR_H1ERROR, derrorH1,&
              rdataSet%Rvector(j), convst_refFunction, rcollection)

          write(*,fmt='(A,T20,G20.10,X,G20.10,X,G20.10)')&
              trim(adjustl(rdataSet%SvarNames(j))), derrorL1, derrorL2, derrorH1

          ! Release temporal 1-block vector
          call lsysbl_releaseVector(rvector)

        end if

      end do
    end do

  end subroutine convst_computeErrorNorms

  !*****************************************************************************

!<subroutine>

  subroutine convst_releaseDataSet(rdataSet)

!<desciption>
    ! This subroutine releases all memory of the given data set structure
!</desciption>

!<inputoutput>
    ! Data set whose memory should be released
    type(t_dataSet), intent(inout) :: rdataSet
!</inputoutput>

    ! local variables
    integer :: i

    call tria_done(rdataSet%rtriangulation)
    call spdiscr_releaseDiscr(rdataSet%rdiscretisation)

    do i = 1, size(rdataSet%Rvector)
      call lsyssc_releaseVector(rdataSet%Rvector(i))
    end do

    deallocate(rdataSet%Rvector)

  end subroutine convst_releaseDataSet

end program convstudy
