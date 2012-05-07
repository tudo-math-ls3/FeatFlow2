!##############################################################################
!# ****************************************************************************
!# <Name> transport_postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all postprocessing routines which are required to
!# solve scalar conservation laws in arbitrary spatial dimensions.
!#
!# The following routines are available:
!#
!# 1.) transp_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 2.) transp_outputStatistics
!#      -> Outputs the application statitics
!#
!# 3.) transp_projectSolution
!#      -> Performs conservative projection of the solution from
!#         a given FE-space to another FE-space
!#
!# </purpose>
!##############################################################################

module transport_postprocessing

#include "../../flagship.h"

!$use omp_lib
  use bilinearformevaluation
  use collection
  use derivatives
  use flagship_basic
  use fsystem
  use genoutput
  use linearformevaluation
  use lineariser
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use scalarpde
  use spatialdiscretisation
  use statistics
  use stdoperators
  use triangulation
  use ucd

  ! Modules from transport model
  use transport_basic
  use transport_callback

  implicit none

  private

  public :: transp_outputSolution
  public :: transp_outputStatistics
  public :: transp_projectSolution

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_outputSolution(rparlist, ssectionName,&
    rproblemLevel, rsolutionPrimal, rsolutionDual, dtime)

!<description>
    ! This subroutine exports the solution vector to file in UCD format
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! OPTIONAL: solution vector for primal problem
    type(t_vectorBlock), intent(in), optional :: rsolutionPrimal

    ! OPTIONAL: solution vector for dual problem
    type(t_vectorBlock), intent(in), optional :: rsolutionDual

    ! OPTIONAL: simulation time
    real(DP), intent(in), optional :: dtime
!</input>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: sucdsolution

    ! persistent variable
    integer, save :: ifilenumber = 1

    ! local variables
    type(t_ucdExport) :: rexport
    type(t_triangulation) :: rtriangulationPrimal,rtriangulationDual
    type(t_blockDiscretisation) :: rdiscretisationPrimal
    type(t_blockDiscretisation) :: rdiscretisationDual
    type(t_vectorBlock) :: rvectorPrimal,rvectorDual
    real(DP), dimension(:,:), allocatable :: DdofCoords
    real(DP), dimension(:), pointer :: p_DdataPrimal,p_DdataDual
    real(DP), dimension(:), pointer :: p_DdofCoords
    integer :: iformatUCD,ilineariseUCD,nrefineUCD,dofCoords,idofe,idim
    logical :: bexportMeshOnly,bdiscontinuous
  
    
    ! Initialisation
    bexportMeshOnly = .true.
    if (present(rsolutionPrimal) .or.&
        present(rsolutionDual)) bexportMeshOnly=.false.

    nullify(p_DdataPrimal, p_DdataDual)

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'dofCoords', dofCoords, 0)
    call parlst_getvalue_string(rparlist, ssectionName,&
                                'output', soutputName)
    call parlst_getvalue_string(rparlist, trim(soutputName),&
                                'sucdsolution', sucdsolution)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'iformatucd', iformatUCD)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'ilineariseucd', ilineariseUCD, UCDEXPORT_STD)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'nrefineucd', nrefineUCD, 0)   

    ! Initialise the UCD exporter
    select case(ilineariseUCD)
    case (UCDEXPORT_STD)
      call flagship_initUCDexport(rproblemLevel,&
          sucdsolution, iformatUCD, rexport, ifilenumber)
      
      ! Set pointers to solution(s)
      if (present(rsolutionPrimal))&
          call lsysbl_getbase_double(rsolutionPrimal, p_DdataPrimal)
      if (present(rsolutionDual))&
          call lsysbl_getbase_double(rsolutionDual, p_DdataDual)
      
    case (UCDEXPORT_P1CONTINUOUS,&
          UCDEXPORT_P1DISCONTINUOUS)
      bdiscontinuous = (ilineariseUCD .eq. UCDEXPORT_P1DISCONTINUOUS)
      
      if (present(rsolutionPrimal)) then
        call lin_lineariseVectorGlobal(rsolutionPrimal, rdiscretisationPrimal,&
            rtriangulationPrimal, rvectorPrimal, nrefineUCD, 0, bdiscontinuous)
        call lsysbl_getbase_double(rvectorPrimal, p_DdataPrimal)
        
        if (present(rsolutionDual)) then
          call lin_lineariseVectorGlobal(rsolutionDual, rdiscretisationDual,&
              rtriangulationDual, rvectorDual, nrefineUCD, 0, bdiscontinuous)
          call lsysbl_getbase_double(rvectorDual, p_DdataDual)
        end if
        
        ! We assume that both primal and dual solutions are based on
        ! the same triangulation, thus both vectors can be exported
        ! using the same triangulation.
        call flagship_initUCDexport(rproblemLevel, sucdsolution,&
            iformatUCD, rexport, ifilenumber, rtriangulationPrimal)

      elseif (present(rsolutionDual)) then
        call lin_lineariseVectorGlobal(rsolutionDual, rdiscretisationDual,&
            rtriangulationDual, rvectorDual, nrefineUCD, 0, bdiscontinuous)
        call lsysbl_getbase_double(rvectorDual, p_DdataDual)
        
        call flagship_initUCDexport(rproblemLevel, sucdsolution,&
            iformatUCD, rexport, ifilenumber, rtriangulationDual)
      end if

    case default
      call output_line('Unsupported type of solution output!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_outputSolution')
      call sys_halt()
    end select

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Prepare array containing the coordinates of the DOFs
    if (.not.bexportMeshOnly .and. dofCoords .gt. 0) then
      ! Get coordinate vector (1D-format)
      call lsyssc_getbase_double(&
          rproblemLevel%RvectorBlock(dofCoords)%RvectorBlock(1), p_DdofCoords)
      
      ! Allocate temporal memory
      allocate(DdofCoords(rproblemLevel%rtriangulation%ndim,&
                          size(p_DdofCoords)/rproblemLevel%rtriangulation%ndim))

      ! Recast coordinates into 2D-format (1:NDIM,1:NDOF)
      do idofe = 1, size(DdofCoords,2)
        do idim = 1, rproblemLevel%rtriangulation%ndim
          DdofCoords(idim,idofe) =&
              p_DdofCoords((idofe-1)*rproblemLevel%rtriangulation%ndim+idim)
        end do
      end do

      ! Set tracer coordinates
      call ucd_setTracers(rexport, DdofCoords)

      ! Release temporal memory
      deallocate(DdofCoords)
    end if
    
    ! Add primal solution vector
    if (associated(p_DdataPrimal)) then
      call ucd_addVariableVertexBased (rexport, 'u',&
          UCD_VAR_STANDARD, p_DdataPrimal)
      
      ! Add pointvalues as tracer variable?
      if (dofCoords .gt. 0)&
          call ucd_addTracerVariable(rexport, 'u', p_DdataPrimal)
    end if

    ! Add dual solution vector
    if (associated(p_DdataDual)) then
      call ucd_addVariableVertexBased (rexport, 'z',&
          UCD_VAR_STANDARD, p_DdataDual)

      ! Add pointvalues as tracer variable?
      if (dofCoords .gt. 0)&
          call ucd_addTracerVariable(rexport, 'z', p_DdataDual)
    end if

    ! Write UCD file
    call ucd_write(rexport)
    call ucd_release(rexport)

    ! Release temporal memory
    call lsysbl_releaseVector(rvectorPrimal)
    call lsysbl_releaseVector(rvectorDual)
    call spdiscr_releaseBlockDiscr(rdiscretisationPrimal)
    call spdiscr_releaseBlockDiscr(rdiscretisationDual)
    call tria_done(rtriangulationPrimal)
    call tria_done(rtriangulationDual)

  end subroutine transp_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine transp_outputStatistics(rtimerTotal, ssectionName, rcollection)

!<description>
    ! This subroutine output application statistics
!</description>

!<input>
    ! timer for total time measurement
    type(t_timer), intent(in) :: rtimerTotal

    ! section name in collection structure
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff
    type(t_timer), pointer :: p_rtimerAssemblyMatrix
    type(t_timer), pointer :: p_rtimerAssemblyVector
    type(t_timer), pointer :: p_rtimerPrePostprocess
    real(DP) :: dtotalTime, dfraction


    ! Get timer objects from collection
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)
    p_rtimerAssemblyMatrix => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', ssectionName=ssectionName)
    p_rtimerAssemblyVector => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)

    ! Output statistics
    call output_lbrk()
    call output_line('Time measurement: '//trim(adjustl(ssectionName)))
    call output_line('-----------------')

    call stat_subTimers(p_rtimerAssemblyMatrix, p_rtimerSolution)
    call stat_subTimers(p_rtimerAssemblyVector, p_rtimerSolution)

    dtotalTime = max(rtimerTotal%delapsedCPU, rtimerTotal%delapsedReal)
    dfraction  = 100.0_DP/dtotalTime

    call output_line('Time for computing solution   : '//&
                     trim(adjustl(sys_sdE(p_rtimerSolution%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerSolution%delapsedCPU, 5)))//' %')
    call output_line('Time for mesh adaptivity      : '//&
                     trim(adjustl(sys_sdE(p_rtimerAdaptation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAdaptation%delapsedCPU, 5)))//' %')
    call output_line('Time for error estimation     : '//&
                     trim(adjustl(sys_sdE(p_rtimerErrorEstimation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerErrorEstimation%delapsedCPU, 5)))//' %')
    call output_line('Time for triangulation        : '//&
                     trim(adjustl(sys_sdE(p_rtimerTriangulation%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerTriangulation%delapsedCPU, 5)))//' %')
    call output_line('Time for coefficient assembly : '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyCoeff%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyCoeff%delapsedCPU, 5)))//' %')
    call output_line('Time for matrix assembly      : '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyMatrix%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyMatrix%delapsedCPU, 5)))//' %')
    call output_line('Time for vector assembly      : '//&
                     trim(adjustl(sys_sdE(p_rtimerAssemblyVector%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerAssemblyVector%delapsedCPU, 5)))//' %')
    call output_line('Time for pre-/post-processing : '//&
                     trim(adjustl(sys_sdE(p_rtimerPrePostprocess%delapsedCPU, 5)))//'  '//&
                     trim(adjustl(sys_sdE(dfraction*p_rtimerPrePostprocess%delapsedCPU, 5)))//' %')
    call output_lbrk()
    call output_line('Time for total simulation     : '//&
                     trim(adjustl(sys_sdE(dtotalTime, 5))))
    call output_lbrk()

  end subroutine transp_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine transp_projectSolution(rsourceVector, rdestVector)
    
!<description>
    ! This subroutine performs conservative projection of the given solution
    ! stored in rsourceVector to another FE-space and stores the result in
    ! rdestVector. An FCT algorithm is used ensure monotonicity preservation.
!</description>

!<input>
    ! Source vector
    type(t_vectorBlock), intent(in), target :: rsourceVector
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: rdestVector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_linearForm) :: rform
    type(t_collection) :: rcollection
    type(t_matrixScalar) :: rmatrix1,rmatrix2
    type(t_vectorScalar) :: rvector
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP) :: dmass1, dmass2
    integer :: iblock
    
    ! Set up the linear form
    rform%itermCount = 1
    rform%Idescriptors(1) = DER_FUNC
    
    ! Set up the collection structure
    call collct_init(rcollection)
    rcollection%IquickAccess(1) = SYSTEM_BLOCKFORMAT
    rcollection%p_rvectorQuickAccess1 => rsourceVector
    
    ! Assemble the linear form for destination vector
    do iblock = 1, rdestVector%nblocks
      
      ! Create sparsity pattern of mass matrices
      call bilf_createMatrixStructure(&
          rsourceVector%p_rblockDiscr%RspatialDiscr(1),&
          LSYSSC_MATRIX9, rmatrix1)
      call bilf_createMatrixStructure(&
          rdestVector%p_rblockDiscr%RspatialDiscr(1),&
          LSYSSC_MATRIX9, rmatrix2)

      ! Create mass matrices
      call stdop_assembleSimpleMatrix(rmatrix1, DER_FUNC, DER_FUNC, 1.0_DP, .true.)
      call stdop_assembleSimpleMatrix(rmatrix2, DER_FUNC, DER_FUNC, 1.0_DP, .true.)
      
      ! Compute the lumped mass matrices
      call lsyssc_lumpMatrixScalar(rmatrix1, LSYSSC_LUMP_DIAG)
      call lsyssc_lumpMatrixScalar(rmatrix2, LSYSSC_LUMP_DIAG)

      ! Set the number of the scalar subvector to the collection structure
      rcollection%IquickAccess(2) = iblock
      
      ! Assemble the linear form for the scalar subvector
      call linf_buildVectorScalar(rform, .true.,&
          rdestVector%RvectorBlock(iblock),&
          fcoeff_buildVectorSc_sim=transp_coeffVectorFE,&
          rcollection=rcollection)

      ! Compute the lumped L2-projection
      call lsyssc_invertedDiagMatVec(rmatrix2, rdestVector%RvectorBlock(iblock),&
          1.0_DP, rdestVector%RvectorBlock(iblock))
      
      ! Compute density-mass
      call lsyssc_duplicateVector(rsourceVector%RvectorBlock(iblock), rvector,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
      call lsyssc_scalarMatVec(rmatrix1, rsourceVector%RvectorBlock(iblock),&
          rvector, 1.0_DP, 0.0_DP)
      call lsyssc_getbase_double(rvector, p_Ddata)
      dmass1 = sum(p_Ddata)
      call lsyssc_releaseVector(rvector)
      
      ! Compute density-mass
      call lsyssc_duplicateVector(rdestVector%RvectorBlock(iblock), rvector,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
      call lsyssc_scalarMatVec(rmatrix2, rdestVector%RvectorBlock(iblock),&
          rvector, 1.0_DP, 0.0_DP)
      call lsyssc_getbase_double(rvector, p_Ddata)
      dmass2 = sum(p_Ddata)
      call lsyssc_releaseVector(rvector)
            
      ! Release matrices
      call lsyssc_releaseMatrix(rmatrix1)
      call lsyssc_releaseMatrix(rmatrix2)

    end do
    
    ! Release the collection structure
    call collct_done(rcollection)

  end subroutine transp_projectSolution

end module transport_postprocessing
