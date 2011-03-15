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

  use collection
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

    ! Create spatial discretisation
    call spdiscr_initDiscr_simple(rdataSet%rdiscretisation, EL_E001_1D,&
        SPDISC_CUB_AUTOMATIC, rdataSet%rtriangulation)

!!$    call spdiscr_initDiscr_triquad(rdataSet%rdiscretisation, EL_E001, EL_E011,&
!!$        SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, rdataSet%rtriangulation)

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

  !*****************************************************************************

!<subroutine>

  subroutine convst_refFunction(cderivative, rdiscretisation, nelements,&
      npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset, Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors. It has to compute
    ! the (analytical) values of a function in a couple of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
!</description>
    
!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
!</input>
  
!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    call convst_evaluateFE(rdiscretisation%ndimension, npointsPerElement*nelements,&
        cderivative, Dvalues, rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints)
    
  end subroutine convst_refFunction

 !*****************************************************************************

!<subroutine>

  subroutine convst_evaluateFE(ndim,npoints, iderType, Dvalues, rvectorScalar,&
      Dpoints, Ielements, IelementsHint, cnonmeshPoints)

!<description>
    ! This is a wrapper routine which allows to call subroutine
    ! fevl_evaluate1 for (1..npointsPerElement,1..nelements) data
!</description>

!<input>
    ! Number of spatial dimension
    integer, intent(in) :: ndim

    ! Number of points
    integer, intent(in) :: npoints

    ! Type of function value to evaluate. One of the DER_xxxx constants,
    ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
    integer, intent(in) :: iderType
    
    ! The scalar solution vector that is to be evaluated.
    type(t_vectorScalar), intent(in) :: rvectorScalar
    
    ! A list of points where to evaluate.
    ! DIMENSION(1..ndim,1..npoints)
    real(DP), dimension(ndim,npoints), intent(in) :: Dpoints

    ! OPTIONAL: A list of elements containing the points Dpoints.
    ! If this is not specified, the element numbers containing the points
    ! are determined automatically.
    integer, dimension(npoints), intent(in), optional :: Ielements

    ! OPTIONAL: A list of elements that are near the points in Dpoints.
    ! This gives only a hint where to start searching for the actual elements
    ! containing the points. This is ignored if Ielements is specified!
    integer, dimension(npoints), intent(in), optional :: IelementsHint
  
    ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens
    ! if a point is located outside of the domain. May happen e.g. in
    ! nonconvex domains. FEVL_NONMESHPTS_NONE is the default 
    ! parameter if cnonmeshPoints is not specified. 
    integer, intent(in), optional :: cnonmeshPoints
!</input>

!<output>
    ! Values of the FE function at the points specified by Dpoints.
    real(DP), dimension(npoints), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! Call FE-evaluation routines
    call fevl_evaluate(iderType, Dvalues, rvectorScalar, Dpoints,&
        Ielements, IelementsHint, cnonmeshPoints)

  end subroutine convst_evaluateFE

end program convstudy
