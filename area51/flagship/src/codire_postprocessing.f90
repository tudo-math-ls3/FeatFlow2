!##############################################################################
!# ****************************************************************************
!# <name> codire_postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to postprocess
!# the solution to a conservation law for a scalar variable
!#
!# The following routines are available:
!#
!# 1.) codire_outputSolution
!#     -> write the solution to file in UCD format
!#
!# 2.) codire_outputVectorScalar
!#     -> write a scalar vector to file in UCD format
!#
!# 3.) codire_outputVectorBlock
!#     -> write a block vector to file in UCD format
!#
!# 4.) codire_calcSolutionError
!#     -> compute the L1- and L2-norm of the solution error
!#        and its standard deviation
!#
!# The following auxiliary routines are available:
!#
!# 1.) codire_getExactSolution
!#     -> compute the values of the exact solution and its derivatives.
!#
!# </purpose>
!##############################################################################

module codire_postprocessing

  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use pprocerror
  use problem
  use statistics
  use storage
  use ucd

  use codire_basic
  use codire_init

  implicit none

  private

  public :: codire_outputSolution
  public :: codire_outputVectorScalar
  public :: codire_outputVectorBlock
  public :: codire_calcSolutionError

contains

  !*****************************************************************************

!<subroutine>

  subroutine codire_outputSolution(rproblemLevel, rsolution, ttime, sfilename,&
                                ioutputUCD, rdualSolution, breset)

!<description>
    ! This subroutine outputs the solution values to UCD file.
!</description>

!<input>
    ! The multigrid structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! The solution vector
    type(t_vectorBlock), intent(IN) :: rsolution

    ! The simulation time
    real(DP), intent(IN) :: ttime

    ! The name of the output file
    character(LEN=*), intent(IN) :: sfilename

    ! The type of UCD output
    integer, intent(IN) :: ioutputUCD

    ! OPTIONAL: the dual solution
    type(t_vectorBlock), intent(IN), optional :: rdualSolution

    ! OPTIONAL: Switch which decides whether new 
    !           enumerator should be started
    logical, intent(IN), optional :: breset
!</input>
!</subroutine>

    ! local variables
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:), pointer :: p_DvelocityX, p_DvelocityY, p_DvelocityZ
    integer :: ienumGMV  = 0
    integer :: ienumAVS  = 0
    integer :: ienumVTK  = 0

    ! Start time measurement
    call stat_startTimer(rtimer_prepostprocess, STAT_TIMERSHORT)

    ! Start UCD export to file
    select case(ioutputUCD)
    case (UCD_FORMAT_GMV,UCD_FORMAT_BGMV)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumGMV = 0
      end if

      ienumGMV = ienumGMV+1
      if (ioutputUCD .eq. UCD_FORMAT_GMV) then
        call ucd_startGMV (rexport, UCD_FLAG_STANDARD,&
                           rproblemLevel%rtriangulation,&
                           trim(adjustl(sfilename))//'.'//&
                           trim(sys_si0(ienumGMV,5))//'.gmv')
      else
        call ucd_startBGMV (rexport, UCD_FLAG_STANDARD,&
                            rproblemLevel%rtriangulation,&
                            trim(adjustl(sfilename))//'.'//&
                            trim(sys_si0(ienumGMV,5))//'.gmv')
      end if

    case (UCD_FORMAT_AVS)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumAVS = 0
      end if

      ienumAVS = ienumAVS+1
      call ucd_startAVS (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//'.'//&
                         trim(sys_si0(ienumAVS,5))//'.avs')

    case (UCD_FORMAT_VTK)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumVTK = 0
      end if

      ienumVTK = ienumVTK+1
      call ucd_startVTK (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//'.'//&
                         trim(sys_si0(ienumVTK,5))//'.vtu')

    case DEFAULT
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_outputSolution')
      call sys_halt()
    end select

    ! Add solution time
    call ucd_setSimulationTime(rexport, ttime)

    ! Add velocity vector
    ! How many spatial dimensions do we have?
    select case(rproblemLevel%rdiscretisation%ndimension)

    case (NDIM1D)
      call lsyssc_getbase_double(&
          rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(1), p_DvelocityX)

      call ucd_addVarVertBasedVec(rexport, 'velocity',&
                                  p_DvelocityX(1:rproblemLevel%rtriangulation%NVT))

    case (NDIM2D)
      call lsyssc_getbase_double(&
          rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(&
          rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(2), p_DvelocityY)
      
      call ucd_addVarVertBasedVec(rexport, 'velocity',&
                                  p_DvelocityX(1:rproblemLevel%rtriangulation%NVT),&
                                  p_DvelocityY(1:rproblemLevel%rtriangulation%NVT))

    case (NDIM3D)
      call lsyssc_getbase_double(&
          rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(&
          rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(2), p_DvelocityY)
      call lsyssc_getbase_double(&
          rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(3), p_DvelocityZ)

      call ucd_addVarVertBasedVec(rexport, 'velocity',&
                                  p_DvelocityX(1:rproblemLevel%rtriangulation%NVT),&
                                  p_DvelocityY(1:rproblemLevel%rtriangulation%NVT),&
                                  p_DvelocityZ(1:rproblemLevel%rtriangulation%NVT))
      
    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                        OU_CLASS_ERROR,OU_MODE_STD,'codire_outputSolution')
      call sys_halt()
    end select

    ! Add solution vector
    call lsysbl_getbase_double(rsolution, p_Ddata)
    call ucd_addVariableVertexBased (rexport, 'u', UCD_VAR_STANDARD, p_Ddata)

    ! Add dual solution vector (if present)
    if (present(rdualSolution)) then
      call lsysbl_getbase_double(rdualSolution, p_Ddata)
      call ucd_addVariableVertexBased (rexport, 'z', UCD_VAR_STANDARD, p_Ddata)
    end if
    
    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

    ! Stop time measurement
    call stat_stopTimer(rtimer_prepostprocess)

  end subroutine codire_outputSolution  

  !*****************************************************************************

!<subroutine>

  subroutine codire_outputVectorScalar(rproblemLevel, rvector, ttime, sfilename,&
                                    ioutputUCD, breset)

!<description>
    ! This subroutine outputs the values of a scalar vector to UCD file.
!</description>

!<input>
    ! The multigrid structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! The scalar vector
    type(t_vectorScalar), intent(IN) :: rvector

    ! The simulation time
    real(DP), intent(IN) :: ttime

    ! The name of the output file
    character(LEN=*), intent(IN) :: sfilename

    ! The type of UCD output
    integer, intent(IN) :: ioutputUCD

    ! OPTIONAL: Switch which decides whether new enumerator should be started
    logical, intent(IN), optional :: breset
!</input>
!</subroutine>

    ! local variables
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: ienumGMV  = 0
    integer :: ienumAVS  = 0
    integer :: ienumVTK  = 0

    ! Start time measurement
    call stat_startTimer(rtimer_prepostprocess, STAT_TIMERSHORT)

    ! Start UCD export to file
    select case(ioutputUCD)
    case (UCD_FORMAT_GMV,UCD_FORMAT_BGMV)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumGMV = 0
      end if

      ienumGMV = ienumGMV+1
      if (ioutputUCD .eq. UCD_FORMAT_GMV) then
        call ucd_startGMV (rexport, UCD_FLAG_STANDARD,&
                           rproblemLevel%rtriangulation,&
                           trim(adjustl(sfilename))//'.gmv'//trim(sys_si0(ienumGMV,5)))
      else
        call ucd_startBGMV (rexport, UCD_FLAG_STANDARD,&
                            rproblemLevel%rtriangulation,&
                            trim(adjustl(sfilename))//'.gmv'//trim(sys_si0(ienumGMV,5)))
      end if

    case (UCD_FORMAT_AVS)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumAVS = 0
      end if

      ienumAVS = ienumAVS+1
      call ucd_startAVS (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//trim(sys_si0(ienumAVS,5))//'.avs')

    case (UCD_FORMAT_VTK)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumVTK = 0
      end if

      ienumVTK = ienumVTK+1
      call ucd_startVTK (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//trim(sys_si0(ienumVTK,5))//'.vtu')

    case DEFAULT
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_outputScalarVector')
      call sys_halt()
    end select

    ! Add solution time
    call ucd_setSimulationTime(rexport, ttime)
    
    ! Set pointer
    call lsyssc_getbase_double(rvector, p_Ddata)

    ! What kind of vector are we
    if (rvector%NEQ .eq. rproblemLevel%rtriangulation%NEL) then
      
      ! Add element based error estimator
      call ucd_addVariableElementBased (rexport, 'elem', UCD_VAR_STANDARD, p_Ddata)

    elseif (rvector%NEQ .eq. rproblemLevel%rtriangulation%NVT) then

      ! Add vertex based error estimator
      call ucd_addVariableVertexBased (rexport, 'vert', UCD_VAR_STANDARD, p_Ddata)

    else
      
      call output_line('Unsupported vector!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_outputScalarVector')
      call sys_halt()
      
    end if
    
    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)
    
    ! Stop time measurement
    call stat_stopTimer(rtimer_prepostprocess)

  end subroutine codire_outputVectorScalar

  !*****************************************************************************

!<subroutine>

  subroutine codire_outputVectorBlock(rproblemLevel, rvector, ttime, sfilename,&
                                   ioutputUCD, breset)

!<description>
    ! This subroutine outputs the values of a block vector to UCD file.
!</description>

!<input>
    ! The multigrid structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! The block vector
    type(t_vectorBlock), intent(IN) :: rvector

    ! The simulation time
    real(DP), intent(IN) :: ttime

    ! The name of the output file
    character(LEN=*), intent(IN) :: sfilename

    ! The type of UCD output
    integer, intent(IN) :: ioutputUCD

    ! OPTIONAL: Switch which decides whether new enumerator should be started
    logical, intent(IN), optional :: breset
!</input>
!</subroutine>

    ! local variables
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: ienumGMV  = 0
    integer :: ienumAVS  = 0
    integer :: ienumVTK  = 0
    integer :: iblock

    ! Start time measurement
    call stat_startTimer(rtimer_prepostprocess, STAT_TIMERSHORT)

    ! Start UCD export to file
    select case(ioutputUCD)
    case (UCD_FORMAT_GMV,UCD_FORMAT_BGMV)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumGMV = 0
      end if

      ienumGMV = ienumGMV+1
      if (ioutputUCD .eq. UCD_FORMAT_GMV) then
        call ucd_startGMV (rexport, UCD_FLAG_STANDARD,&
                           rproblemLevel%rtriangulation,&
                           trim(adjustl(sfilename))//'.gmv'//trim(sys_si0(ienumGMV,5)))
      else
        call ucd_startBGMV (rexport, UCD_FLAG_STANDARD,&
                            rproblemLevel%rtriangulation,&
                            trim(adjustl(sfilename))//'.gmv'//trim(sys_si0(ienumGMV,5)))
      end if

    case (UCD_FORMAT_AVS)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumAVS = 0
      end if

      ienumAVS = ienumAVS+1
      call ucd_startAVS (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//trim(sys_si0(ienumAVS,5))//'.avs')

    case (UCD_FORMAT_VTK)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumVTK = 0
      end if

      ienumVTK = ienumVTK+1
      call ucd_startVTK (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//trim(sys_si0(ienumVTK,5))//'.vtu')

    case DEFAULT
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_outputBlockVector')
      call sys_halt()
    end select

    ! Add solution time
    call ucd_setSimulationTime(rexport, ttime)
    
    ! Loop over all blocks
    do iblock = 1, rvector%nblocks
      
      ! Set pointer
      call lsyssc_getbase_double(rvector%Rvectorblock(iblock), p_Ddata)
      
      ! What kind of vector are we
      if (rvector%RvectorBlock(iblock)%NEQ .eq. rproblemLevel%rtriangulation%NEL) then
        
        ! Add element based error estimator
        call ucd_addVariableElementBased (rexport, 'elem'//trim(sys_si(iblock,2)),&
                                          UCD_VAR_STANDARD, p_Ddata)
        
      elseif (rvector%NEQ .eq. rproblemLevel%rtriangulation%NVT) then
        
        ! Add vertex based error estimator
        call ucd_addVariableVertexBased (rexport, 'vert'//trim(sys_si(iblock,2)),&
                                         UCD_VAR_STANDARD, p_Ddata)
        
      else
        
        call output_line('Unsupported vector!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_outputBlockVector')
        call sys_halt()
        
      end if

    end do
    
    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)
    
    ! Stop time measurement
    call stat_stopTimer(rtimer_prepostprocess)

  end subroutine codire_outputVectorBlock

  !*****************************************************************************

!<subroutine>

  subroutine codire_calcSolutionError(rproblemLevel, rsolution, indatfile,&
                                   ttime, rcollection)

!<description>
    ! This subroutine computes the $L_1$- and $L_2$-error of the
    ! solution error $u-u_h$. If no exact solution is available, then
    ! this subroutine returns immediately.
!</description>

!<input>
    ! multigrid level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution

    ! name if ondat file
    character(LEN=*), intent(IN) :: indatfile

    ! simulation time
    real(DP), intent(IN) :: ttime
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure. This structure is given to the
    ! callback function to provide additional information.)
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local constants
    character(LEN=*), dimension(NDIM1D+1), parameter ::&
                      SOL_SYMBOLICVARS1D = (/ (/'x'/), (/'t'/) /)
    character(LEN=*), dimension(NDIM2D+1), parameter ::&
                      SOL_SYMBOLICVARS2D = (/ (/'x'/), (/'y'/), (/'t'/) /)
    character(LEN=*), dimension(NDIM3D+1), parameter ::&
                      SOL_SYMBOLICVARS3D = (/ (/'x'/), (/'y'/), (/'z'/), (/'t'/) /)

    ! local variables
    type(t_problem), pointer :: p_rproblem
    type(t_fparser) :: rparser
    type(t_vectorBlock) :: rsolutionExact
    real(DP) :: derrorL1,derrorL2,ddeviation,ddeviationExact
    integer :: istatus
    
    
    ! Initialize exact solution
    call codire_initExactSolution(rproblemLevel, rsolutionExact, indatfile, ttime, istatus)

    ! If no exact solution is available, return
    if (istatus .ne. 0) return

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem

    ! How many spatial dimensions do we have?
    select case(rproblemLevel%rdiscretisation%ndimension)
      
    case (NDIM1D)
      ! Create solution profile in 1D
      call problem_createProfile(p_rproblem, indatfile, '[exact_solution]',&
                                 SOL_SYMBOLICVARS1D, rparser, istatus, ttime)

    case (NDIM2D)
      ! Create solution profile in 2D
      call problem_createProfile(p_rproblem, indatfile, '[exact_solution]',&
                                 SOL_SYMBOLICVARS2D, rparser, istatus, ttime)
      
    case (NDIM3D)
      ! Create solution profile in 3D
      call problem_createProfile(p_rproblem, indatfile, '[exact_solution]',&
                                 SOL_SYMBOLICVARS3D, rparser, istatus, ttime)
      
    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcSolutionError')
      call sys_halt()
    end select

    ! Add function parser to collection so that it 
    ! can be used from within the callback functions
    call collct_setvalue_pars(rcollection, 'parser_exactsolution', rparser, .true.)

    ! Add time to quick access field
    rcollection%Dquickaccess(1) = ttime

    ! Add dimension to quick access field
    rcollection%Iquickaccess(1) = rproblemLevel%rdiscretisation%ndimension


    ! Compute L1-errors
    call pperr_scalarErrorEstimate(rsolution%RvectorBlock(1),&
                                   rsolutionExact%RvectorBlock(1),&
                                   PPERR_L1ERROR, derrorL1)
    call output_line('L1-error:                 '//trim(sys_sdE(derrorL1,16)))

    call pperr_scalar(rsolution%RvectorBlock(1), PPERR_L1ERROR, derrorL1,&
        codire_getExactSolution, rcollection)
    call output_line('L1-error:                 '//trim(sys_sdE(derrorL1,16)))


    ! Compute L2-errors
    call pperr_scalarErrorEstimate(rsolution%RvectorBlock(1),&
                                   rsolutionExact%RvectorBlock(1),&
                                   PPERR_L2ERROR, derrorL2)
    call output_line('L2-error:                 '//trim(sys_sdE(derrorL2,16)))

    call pperr_scalar(rsolution%RvectorBlock(1), PPERR_L2ERROR, derrorL2,&
                      codire_getExactSolution, rcollection)
    call output_line('L2-error:                 '//trim(sys_sdE(derrorL2,16)))


    ! Compute standard deviation
    call pperr_scalarStandardDeviation(rsolution%RvectorBlock(1),&
                                       ddeviation)
    call output_line('Standard deviation:       '//trim(sys_sdE(ddeviation,16)))

    
    ! Compute deviation of exact solution
    call pperr_scalarStandardDeviation(rsolutionExact%RvectorBlock(1),&
                                       ddeviationExact)
    call output_line('Relative deviation error: '//&
        trim(sys_sdE((ddeviation**2-ddeviationExact**2)/&
        max(SYS_EPSREAL,ddeviationExact**2),16)))
    
    ! Delete function parser from collection
    call collct_deletevalue(rcollection, 'parser_exactsolution')

    ! Release memory
    call fparser_release(rparser)
    call lsysbl_releaseVector(rsolutionExact)
  end subroutine codire_calcSolutionError

  ! ***************************************************************************

!<subroutine>

  subroutine codire_getExactSolution(cderivative, rdiscretisation, nelements,&
                                  npointsPerElement, Dpoints, IdofsTest,&
                                  rdomainIntSubset, Dvalues, rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
  
!<description>
    ! This subroutine is called during the calculation of errors. It has to compute
    ! the (analytical) values of a function in a couple of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
!</description>
  
!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(IN) :: cderivative
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
    
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional :: rcollection
  
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: rparser
    real(DP), dimension(NDIM1D+1) :: DvariableValues1D
    real(DP), dimension(NDIM2D+1) :: DvariableValues2D
    real(DP), dimension(NDIM3D+1) :: DvariableValues3D
    real(DP) :: time
    integer :: i,j,ndimension

    ! Get function parser from collection
    rparser => collct_getvalue_pars(rcollection, 'parser_exactsolution')

    ! Get time and dimension from quick access fields
    time       = rcollection%Dquickaccess(1)
    ndimension = rcollection%Iquickaccess(1)
    
    select case (ndimension)

    case (NDIM1D)

      select case (cderivative)
      case (DER_FUNC)
        
        do j = 1, size(Dvalues, 2)
          do i = 1, size(Dvalues, 1)
            
            DvariableValues1D(1) = Dpoints(1,i,j)
            DvariableValues1D(2) = time
            
            call fparser_evalFunction(rparser, 1, DvariableValues1D, Dvalues(i,j))
          end do
        end do

      case DEFAULT
        ! Unknown. Set the result to 0.0.
        Dvalues = 0.0_DP
      end select

    case (NDIM2D)

      select case (cderivative)
      case (DER_FUNC)
        
        do j = 1, size(Dvalues, 2)
          do i = 1, size(Dvalues, 1)
            
            DvariableValues2D(1) = Dpoints(1,i,j)
            DvariableValues2D(2) = Dpoints(2,i,j)
            DvariableValues2D(3) = time
            
            call fparser_evalFunction(rparser, 1, DvariableValues2D, Dvalues(i,j))
          end do
        end do

      case DEFAULT
        ! Unknown. Set the result to 0.0.
        Dvalues = 0.0_DP
      end select
      
    case (NDIM3D)

      select case (cderivative)
      case (DER_FUNC)
        
        do j = 1, size(Dvalues, 2)
          do i = 1, size(Dvalues, 1)
            
            DvariableValues3D(1) = Dpoints(1,i,j)
            DvariableValues3D(2) = Dpoints(2,i,j)
            DvariableValues3D(3) = Dpoints(3,i,j)
            DvariableValues3D(4) = time
            
            call fparser_evalFunction(rparser, 1, DvariableValues3D, Dvalues(i,j))
          end do
        end do

      case DEFAULT
        ! Unknown. Set the result to 0.0.
        Dvalues = 0.0_DP
      end select

    case DEFAULT
      ! Unknown. Set the result to 0.0.
      Dvalues = 0.0_DP
    end select

  end subroutine codire_getExactSolution
end module codire_postprocessing
