!##############################################################################
!# ****************************************************************************
!# <name> flagship_postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to postprocess
!# the solution to the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!# 1.) fs_outputSolution
!#     -> write the solution to file in UCD format
!#
!# 2.) fs_outputVectorScalar
!#     -> write a scalar vector to file in UCD format
!#
!# 3.) fs_outputVectorBlock
!#     -> write a block vector to file in UCD format
!#
!# 4.) fs_calcSolutionError
!#     -> compute the L1- and L2-norm of the solution error
!#        and its standard deviation
!#
!# </purpose>
!##############################################################################

module flagship_postprocessing

  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use pprocerror
  use statistics
  use storage
  use ucd

  use flagship_basic
  use flagship_init
  use problem

  implicit none

  private

  public :: fs_outputSolution
  public :: fs_outputVectorScalar
  public :: fs_outputVectorBlock
  public :: fs_calcSolutionError

contains

  !*****************************************************************************

!<subroutine>

  subroutine fs_outputSolution(rproblemLevel, rsolution, ttime, sfilename,&
                               ioutputUCD, breset)

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

    ! OPTIONAL: Switch which decides whether new enumerator should be started
    logical, intent(IN), optional :: breset
!</input>
!</subroutine>

    ! local variables
    type(t_ucdExport) :: rexport
    type(t_vectorScalar) :: rvector1,rvector2,rvector3
    real(DP), dimension(:), pointer :: p_Dsolution
    real(DP), dimension(:), pointer :: p_Ddata,p_DvelocityX,p_DvelocityY,p_DvelocityZ
    integer :: ienumGMV = 0
    integer :: ienumAVS = 0
    integer :: ienumVTK = 0

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
                           trim(adjustl(sfilename))//'.gmv'//sys_si0(ienumGMV,5))
      else
        call ucd_startBGMV (rexport, UCD_FLAG_STANDARD,&
                            rproblemLevel%rtriangulation,&
                            trim(adjustl(sfilename))//'.gmv'//sys_si0(ienumGMV,5))
      end if

    case (UCD_FORMAT_AVS)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumAVS = 0
      end if

      ienumAVS = ienumAVS+1
      call ucd_startAVS (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//sys_si0(ienumAVS,5)//'.avs')

    case (UCD_FORMAT_VTK)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumVTK = 0
      end if

      ienumVTK = ienumVTK+1
      call ucd_startVTK (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//sys_si0(ienumVTK,5)//'.vtu')

    case DEFAULT
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_outputSolution')
      call sys_halt()
    end select

    ! Add solution time
    call ucd_setSimulationTime(rexport, ttime)
    
    ! Set pointers
    call lsysbl_getbase_double(rsolution, p_Dsolution)
    
    ! Create auxiliary vectors
    select case(rproblemLevel%rdiscretisation%ndimension)
    case (NDIM1D)
      call lsyssc_createVector(rvector1, int(size(p_Dsolution)/NVAR1D, I32), .false.)
      call lsyssc_getbase_double(rvector1, p_DvelocityX)

    case (NDIM2D)
      call lsyssc_createVector(rvector1, int(size(p_Dsolution)/NVAR2D, I32), .false.)
      call lsyssc_createVector(rvector2, int(size(p_Dsolution)/NVAR2D, I32), .false.)
      call lsyssc_getbase_double(rvector1, p_DvelocityX)
      call lsyssc_getbase_double(rvector2, p_DvelocityY)

    case (NDIM3D)
      call lsyssc_createVector(rvector1, int(size(p_Dsolution)/NVAR3D, I32), .false.)
      call lsyssc_createVector(rvector2, int(size(p_Dsolution)/NVAR3D, I32), .false.)
      call lsyssc_createVector(rvector3, int(size(p_Dsolution)/NVAR3D, I32), .false.)
      call lsyssc_getbase_double(rvector1, p_DvelocityX)
      call lsyssc_getbase_double(rvector2, p_DvelocityY)
      call lsyssc_getbase_double(rvector3, p_DvelocityZ)

    case DEFAULT
      call output_line('Invalid number of spatial dimensions',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_outputSolution')
      call sys_halt()
    end select
    call lsyssc_getbase_double(rvector1, p_Ddata)
    
    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      select case(rproblemLevel%rdiscretisation%ndimension)
      case (NDIM1D)
        ! Compute velocity in 1D
        call fs_getVariableNodewise(rvector1%NEQ, NVAR1D, p_Dsolution,&
                                    VAR_MOMENTUM_X, p_DvelocityX)
        
        ! Add velocity vectors in 1D
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_DvelocityX)

        ! Compute additional variables in 1D
        call fs_getVariableNodewise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_DENSITY, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata)
      
        call fs_getVariableNodewise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_MOMENTUM_X+1, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'tot_energy', UCD_VAR_STANDARD, p_Ddata)
      
        call fs_getVariableNodewise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_MOMENTUM_X+2, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata)

        call fs_getVariableNodewise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_MOMENTUM_X+3, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'Mach_numb', UCD_VAR_STANDARD, p_Ddata)


      case (NDIM2D)
        ! Compute velocity in 2D
        call fs_getVariableNodewise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_X, p_DvelocityX)
        call fs_getVariableNodewise(rvector2%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_Y, p_DvelocityY)
        
        ! Add velocity vectors in 2D
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_DvelocityX, p_DvelocityY)

        ! Compute additional variables in 2D
        call fs_getVariableNodewise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_DENSITY, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata)
      
        call fs_getVariableNodewise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_Y+1, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'tot_energy', UCD_VAR_STANDARD, p_Ddata)
      
        call fs_getVariableNodewise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_Y+2, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata)

        call fs_getVariableNodewise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_Y+3, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'Mach_numb', UCD_VAR_STANDARD, p_Ddata)
        

      case (NDIM3D)
        ! Compute velocity in 3D
        call fs_getVariableNodewise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_X, p_DvelocityX)
        call fs_getVariableNodewise(rvector2%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Y, p_DvelocityY)
        call fs_getVariableNodewise(rvector3%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Z, p_DvelocityZ)
        
        ! Add velocity vectors in 3D
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_DvelocityX, p_DvelocityY, p_DvelocityZ)

        ! Compute additional variables in 3D
        call fs_getVariableNodewise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_DENSITY, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata)
      
        call fs_getVariableNodewise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Z+1, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'tot_energy', UCD_VAR_STANDARD, p_Ddata)
      
        call fs_getVariableNodewise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Z+2, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata)

        call fs_getVariableNodewise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Z+3, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'Mach_numb', UCD_VAR_STANDARD, p_Ddata)
      end select


    case (SYSTEM_BLOCKFORMAT)

      select case(rproblemLevel%rdiscretisation%ndimension)
      case (NDIM1D)
        ! Compute velocity in 1D
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR1D, p_Dsolution,&
                                     VAR_MOMENTUM_X, p_DvelocityX)
        
        ! Add velocity vectors in 1D
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_DvelocityX)

        ! Compute additional variables in 1D
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_DENSITY, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_MOMENTUM_X+1, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'tot_energy', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_MOMENTUM_X+2, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR1D, p_Dsolution, VAR_MOMENTUM_X+3, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'Mach_numb', UCD_VAR_STANDARD, p_Ddata)


      case (NDIM2D)
        ! Compute velocity in 2D
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR2D, p_Dsolution,&
                                     VAR_MOMENTUM_X, p_DvelocityX)
        call fs_getVariableBlockwise(rvector2%NEQ, NVAR2D, p_Dsolution,&
                                     VAR_MOMENTUM_Y, p_DvelocityY)
        
        ! Add velocity vectors in 2D
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_DvelocityX, p_DvelocityY)

        ! Compute additional variables in 2D
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_DENSITY, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_Y+1, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'tot_energy', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_Y+2, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR2D, p_Dsolution, VAR_MOMENTUM_Y+3, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'Mach_numb', UCD_VAR_STANDARD, p_Ddata)
        

      case (NDIM3D)
        ! Compute velocity in 3D
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR3D, p_Dsolution,&
                                     VAR_MOMENTUM_X, p_DvelocityX)
        call fs_getVariableBlockwise(rvector2%NEQ, NVAR3D, p_Dsolution,&
                                     VAR_MOMENTUM_Y, p_DvelocityY)
        call fs_getVariableBlockwise(rvector3%NEQ, NVAR3D, p_Dsolution,&
                                     VAR_MOMENTUM_Z, p_DvelocityZ)
        
        ! Add velocity vectors in 3D
        call ucd_addVarVertBasedVec(rexport, 'velocity', p_DvelocityX, p_DvelocityY, p_DvelocityZ)

        ! Compute additional variables in 3D
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_DENSITY, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'density', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Z+1, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'tot_energy', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Z+2, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'pressure', UCD_VAR_STANDARD, p_Ddata)
        
        call fs_getVariableBlockwise(rvector1%NEQ, NVAR3D, p_Dsolution, VAR_MOMENTUM_Z+3, p_Ddata)
        call ucd_addVariableVertexBased (rexport, 'Mach_numb', UCD_VAR_STANDARD, p_Ddata)
      end select


    case DEFAULT
      call output_line('Unsupported system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_outputSolution')
      call sys_halt()
    end select
    
    ! Release temporal memory
    call lsyssc_releaseVector(rvector1)
    call lsyssc_releaseVector(rvector2)
    call lsyssc_releaseVector(rvector3)

    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)

    ! Stop time measurement
    call stat_stopTimer(rtimer_prepostprocess)

  end subroutine fs_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine fs_outputVectorScalar(rproblemLevel, rvector, ttime, sfilename,&
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
    integer :: ienumGMV = 0
    integer :: ienumAVS = 0
    integer :: ienumVTK = 0

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
                           trim(adjustl(sfilename))//'.gmv'//sys_si0(ienumGMV,5))
      else
        call ucd_startBGMV (rexport, UCD_FLAG_STANDARD,&
                            rproblemLevel%rtriangulation,&
                            trim(adjustl(sfilename))//'.gmv'//sys_si0(ienumGMV,5))
      end if

    case (UCD_FORMAT_AVS)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumAVS = 0
      end if

      ienumAVS = ienumAVS+1
      call ucd_startAVS (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//sys_si0(ienumAVS,5)//'.avs')

    case (UCD_FORMAT_VTK)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumVTK = 0
      end if

      ienumVTK = ienumVTK+1
      call ucd_startVTK (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//sys_si0(ienumVTK,5)//'.vtu')

    case DEFAULT
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_outputScalarVector')
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
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_outputScalarVector')
      call sys_halt()
      
    end if
    
    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)
    
    ! Stop time measurement
    call stat_stopTimer(rtimer_prepostprocess)

  end subroutine fs_outputVectorScalar

  !*****************************************************************************

!<subroutine>

  subroutine fs_outputVectorBlock(rproblemLevel, rvector, ttime, sfilename,&
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
    integer :: ienumGMV = 0
    integer :: ienumAVS = 0
    integer :: ienumVTK = 0
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
                           trim(adjustl(sfilename))//'.gmv'//sys_si0(ienumGMV,5))
      else
        call ucd_startBGMV (rexport, UCD_FLAG_STANDARD,&
                            rproblemLevel%rtriangulation,&
                            trim(adjustl(sfilename))//'.gmv'//sys_si0(ienumGMV,5))
      end if

    case (UCD_FORMAT_AVS)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumAVS = 0
      end if

      ienumAVS = ienumAVS+1
      call ucd_startAVS (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//sys_si0(ienumAVS,5)//'.avs')

    case (UCD_FORMAT_VTK)
      ! Do we have to reset the enumerator?
      if (present(breset)) then
        if (breset) ienumVTK = 0
      end if

      ienumVTK = ienumVTK+1
      call ucd_startVTK (rexport, UCD_FLAG_STANDARD,&
                         rproblemLevel%rtriangulation,&
                         trim(adjustl(sfilename))//sys_si0(ienumVTK,5)//'.vtu')

    case DEFAULT
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_outputBlockVector')
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
                         OU_CLASS_ERROR,OU_MODE_STD,'fs_outputBlockVector')
        call sys_halt()
        
      end if

    end do
    
    ! Write UCD file
    call ucd_write  (rexport)
    call ucd_release(rexport)
    
    ! Stop time measurement
    call stat_stopTimer(rtimer_prepostprocess)

  end subroutine fs_outputVectorBlock

  !*****************************************************************************

!<subroutine>

  subroutine fs_calcSolutionError(rproblemLevel, rsolution, indatfile,&
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

    ! local variables
    type(t_problem), pointer :: p_rproblem
    type(t_fparser) :: rparser
    type(t_vectorBlock) :: rsolutionExact
    type(t_vectorScalar) :: rsolutionScalar,rsolutionExactScalar
    real(DP) :: derrorL1,derrorL2,ddeviation,ddeviationExact
    integer :: istatus,ivar
    
    
    ! Initialize exact solution
    call fs_initExactSolution(rproblemLevel, rsolutionExact,&
                              indatfile, ttime, istatus)

    ! If no exact solution is available, return
    if (istatus .ne. 0) return

    ! Set pointer
    p_rproblem => rproblemLevel%p_rproblem
    
    ! What system format are we?
    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      ! Create temporal vectors for (exact) solution
      call lsyssc_createVecByDiscr(&
          rsolution%p_rblockDiscr%RspatialDiscr(1),&
          rsolutionScalar, .false., rsolution%cdataType)
      
      call lsyssc_createVecByDiscr(&
          rsolutionExact%p_rblockDiscr%RspatialDiscr(1),&
          rsolutionExactScalar, .false., rsolutionExact%cdataType)

      do ivar = 1, fs_getNVAR(rproblemLevel)
        call output_line('Variable '//trim(sys_si(ivar,2)))
        
        call lsyssc_packVector(rsolution%RvectorBlock(1),&
                               rsolutionScalar, ivar)
        call lsyssc_packVector(rsolutionExact%RvectorBlock(1),&
                               rsolutionExactScalar, ivar)

        ! Compute L1-errors
        call pperr_scalarErrorEstimate(rsolutionScalar, rsolutionExactScalar,&
                                       PPERR_L1ERROR, derrorL1)
        call output_line('L1-error:                 '//trim(sys_sdE(derrorL1,16)))
        
        ! Compute L2-errors
        call pperr_scalarErrorEstimate(rsolutionScalar, rsolutionExactScalar,&
                                       PPERR_L2ERROR, derrorL2)
        call output_line('L2-error:                 '//trim(sys_sdE(derrorL2,16)))

        ! Compute standard deviation
        call pperr_scalarStandardDeviation(rsolutionScalar, ddeviation)
        call output_line('Standard deviation:       '//trim(sys_sdE(ddeviation,16)))

        ! Compute deviation of exact solution
        call pperr_scalarStandardDeviation(rsolutionExactScalar, ddeviationExact)
        call output_line('Relative deviation error: '//&
            trim(sys_sdE((ddeviation**2-ddeviationExact**2)/&
            max(SYS_EPSREAL,ddeviationExact**2),16)))
      end do
            
      ! Release memory
      call lsyssc_releaseVector(rsolutionScalar)
      call lsyssc_releaseVector(rsolutionExactScalar)

    case (SYSTEM_BLOCKFORMAT)

      do ivar = 1, fs_getNVAR(rproblemLevel)
        call output_line('Variable '//trim(sys_si(ivar,2)))

        ! Compute L1-errors
        call pperr_scalarErrorEstimate(rsolution%RvectorBlock(ivar),&
            rsolutionExact%RvectorBlock(ivar), PPERR_L1ERROR, derrorL1)
        call output_line('L1-error:                 '//trim(sys_sdE(derrorL1,16)))
        
        ! Compute L2-errors
        call pperr_scalarErrorEstimate(rsolution%RvectorBlock(ivar),&
        rsolutionExact%RvectorBlock(ivar), PPERR_L2ERROR, derrorL2)
        call output_line('L2-error:                 '//trim(sys_sdE(derrorL2,16)))

        ! Compute standard deviation
        call pperr_scalarStandardDeviation(rsolution%RvectorBlock(ivar),&
            ddeviation)
        call output_line('Standard deviation:       '//trim(sys_sdE(ddeviation,16)))

        ! Compute deviation of exact solution
        call pperr_scalarStandardDeviation(rsolutionExact%RvectorBlock(ivar),&
            ddeviationExact)
        call output_line('Relative deviation error: '//&
            trim(sys_sdE((ddeviation**2-ddeviationExact**2)/&
            max(SYS_EPSREAL,ddeviationExact**2),16)))
      end do

    case DEFAULT
      call output_line('Unsupported system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fs_calcSolutionError')
      call sys_halt()
    end select
    
    ! Release memory
    call fparser_release(rparser)
    call lsysbl_releaseVector(rsolutionExact)
  end subroutine fs_calcSolutionError

end module flagship_postprocessing
