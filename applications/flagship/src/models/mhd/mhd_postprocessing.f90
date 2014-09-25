!##############################################################################
!# ****************************************************************************
!# <Name> mhd_postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all postprocessing routines which are required to
!# solve the compressible ideal MHD equations in arbitrary spatial dimensions
!#
!# The following routines are available:
!#
!# 1.) mhd_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 2.) mhd_outputStatistics
!#      -> Outputs the application statitics
!#
!# 3.) mhd_projectSolution
!#      -> Performs conservative projection of the solution from
!#         a given FE-space to another FE-space
!#
!# 4.) mhd_outputFeSolution
!#     -> Outputs the FE-solution in raw format
!#
!# </purpose>
!##############################################################################

module mhd_postprocessing

#include "flagship.h"
#include "mhd.h"

!$ use omp_lib
  use basicgeometry
  use bilinearformevaluation
  use collection
  use derivatives
  use flagship_basic
  use fsystem
  use genoutput
  use io
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
  use vectorio

  ! Modules from MHD model
  use mhd_basic
  use mhd_basic1d
  use mhd_basic2d
  use mhd_basic3d
  use mhd_callback

  implicit none

  private

  public :: mhd_outputSolution
  public :: mhd_outputStatistics
  public :: mhd_projectSolution
  public :: mhd_outputFeSolution

contains

  !*****************************************************************************

!<subroutine>

  subroutine mhd_outputSolution(rparlist, ssectionName,&
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
    real(DP), dimension(1) :: DtimeAux
    real(DP) :: density_ref, velocity_ref, length_ref
    integer :: isystemFormat,iformatUCD,ilineariseUCD,nrefineUCD
    integer :: dofCoords,idofe,idim
    logical :: bexportMeshOnly,bdiscontinuous
    logical :: bconvert

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
                             'iformatucd', iformatUCD, 0)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'isystemformat', isystemformat)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'ilineariseucd', ilineariseUCD, UCDEXPORT_STD)
    call parlst_getvalue_int(rparlist, trim(soutputName),&
                             'nrefineucd', nrefineUCD, 0)

    ! Get reference values for density, velocity and length to convert
    ! non-dimensional quantities into physical quantities
    call parlst_getvalue_double(rparlist, trim(soutputName),&
                                'density_ref', density_ref, 1.0_DP)
    call parlst_getvalue_double(rparlist, trim(soutputName),&
                                'velocity_ref', velocity_ref, 1.0_DP)
    call parlst_getvalue_double(rparlist, trim(soutputName),&
                                'length_ref', length_ref, 1.0_DP)

    ! Do we have to convert into physical quantities?
    bconvert = ((density_ref  .ne. 1.0_DP) .or.&
                (velocity_ref .ne. 1.0_DP) .or.&
                (length_ref   .ne. 1.0_DP))

    if (iformatUCD .eq. 0) then
      call output_line('No valid output format is specified!',&
          OU_CLASS_WARNING,OU_MODE_STD,'mhd_outputSolution')
      return
    end if

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
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_outputSolution')
      call sys_halt()
    end select

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) then
      DtimeAux = dtime
      if (bconvert) call mhd_convertVariableDim(DtimeAux, 'time',&
          density_ref, velocity_ref, length_ref)
      call ucd_setSimulationTime(rexport, DtimeAux(1))
    end if

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
    if (associated(p_DdataPrimal))&
        call outputSolution(rexport, p_DdataPrimal,'', (dofCoords.gt.0), bconvert)
    
    ! Add dual solution vector
    if (associated(p_DdataDual))&
        call outputSolution(rexport, p_DdataDual,'_dual', (dofCoords.gt.0), bconvert)

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

  contains

    ! Here, the working routine follows

    !***************************************************************************
    ! This subroutine outputs the solution given by the array Ddata

    subroutine outputSolution(rexport, Ddata, csuffix, btracers, bconvert)
      
      ! Input parameters
      real(DP), dimension(:), intent(in) :: Ddata
      character(len=*), intent(in) :: csuffix
      logical, intent(in) :: btracers,bconvert

      ! Input/output paramters
      type(t_ucdExport), intent(inout) :: rexport

      ! local variables
      type(t_vectorScalar) :: rvector1,rvector2,rvector3
      real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
      character(len=SYS_NAMELEN) :: sucdvariable
      integer :: isize,ndim,ivariable,nvariable

      ! Set pointers
      isize = size(Ddata)/mhd_getNVAR(rproblemLevel)
      ndim  = rproblemLevel%rtriangulation%ndim

      ! Create auxiliary vectors
      call lsyssc_createVector(rvector1, isize, .false.)
      call lsyssc_createVector(rvector2, isize, .false.)
      call lsyssc_createVector(rvector3, isize, .false.)
      call lsyssc_getbase_double(rvector1, p_Ddata1)
      call lsyssc_getbase_double(rvector2, p_Ddata2)
      call lsyssc_getbase_double(rvector3, p_Ddata3)

      ! Get number of variables to be written
      nvariable = max(1,&
          parlst_querysubstrings(rparlist,&
          trim(soutputName), 'sucdvariable'))
      
      select case(isystemFormat)
      case(SYSTEM_INTERLEAVEFORMAT)
        
        ! Loop over all variables
        do ivariable = 1, nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sucdvariable', sucdvariable, isubstring=ivariable)
          
          if (trim(sucdvariable) .eq. 'velocity') then
            
            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarInterleaveFormat1d(rvector1%NEQ, NVAR1D,&
                  'velocity_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat1d(rvector2%NEQ, NVAR1D,&
                  'velocity_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat1d(rvector3%NEQ, NVAR1D,&
                  'velocity_z', Ddata, p_Ddata3)
              
            case (NDIM2D)
              call mhd_getVarInterleaveFormat2d(rvector1%NEQ, NVAR2D,&
                  'velocity_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat2d(rvector2%NEQ, NVAR2D,&
                  'velocity_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat2d(rvector3%NEQ, NVAR2D,&
                  'velocity_z', Ddata, p_Ddata3)

            case (NDIM3D)
              call mhd_getVarInterleaveFormat3d(rvector1%NEQ, NVAR3D,&
                  'velocity_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat3d(rvector2%NEQ, NVAR3D,&
                  'velocity_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat3d(rvector3%NEQ, NVAR3D,&
                  'velocity_z', Ddata, p_Ddata3)
            end select

            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, 'velocity',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata2, 'velocity',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata3, 'velocity',&
                  density_ref, velocity_ref, length_ref)
            end if
            
            call ucd_addVarVertBasedVec(rexport, 'velocity'//csuffix,&
                UCD_VAR_VELOCITY, p_Ddata1, p_Ddata2, p_Ddata3)
            
            if (btracers) then
              call ucd_addTracerVariable(rexport, 'velocity_x'//csuffix,&
                  p_Ddata1)
              call ucd_addTracerVariable(rexport, 'velocity_y'//csuffix,&
                  p_Ddata2)
              call ucd_addTracerVariable(rexport, 'velocity_z'//csuffix,&
                  p_Ddata3)
            end if

          elseif (trim(sucdvariable) .eq. 'momentum') then
            
            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarInterleaveFormat1d(rvector1%NEQ, NVAR1D,&
                  'momentum_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat1d(rvector2%NEQ, NVAR1D,&
                  'momentum_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat1d(rvector3%NEQ, NVAR1D,&
                  'momentum_z', Ddata, p_Ddata3)
              
            case (NDIM2D)
              call mhd_getVarInterleaveFormat2d(rvector1%NEQ, NVAR2D,&
                  'momentum_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat2d(rvector2%NEQ, NVAR2D,&
                  'momentum_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat2d(rvector3%NEQ, NVAR2D,&
                  'momentum_z', Ddata, p_Ddata3)

            case (NDIM3D)
              call mhd_getVarInterleaveFormat3d(rvector1%NEQ, NVAR3D,&
                  'momentum_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat3d(rvector2%NEQ, NVAR3D,&
                  'momentum_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat3d(rvector3%NEQ, NVAR3D,&
                  'momentum_z', Ddata, p_Ddata3)
            end select
           
            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, 'momentum',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata2, 'momentum',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata3, 'momentum',&
                  density_ref, velocity_ref, length_ref)
            end if
	    
            call ucd_addVarVertBasedVec(rexport, 'momentum'//csuffix,&
                p_Ddata1, p_Ddata2, p_Ddata3)

            if (btracers) then
              call ucd_addTracerVariable(rexport, 'momentum_x'//csuffix,&
                  p_Ddata1)
              call ucd_addTracerVariable(rexport, 'momentum_y'//csuffix,&
                  p_Ddata2)
              call ucd_addTracerVariable(rexport, 'momentum_z'//csuffix,&
                  p_Ddata3)
            end if
            
          elseif (trim(sucdvariable) .eq. 'magneticfield') then

            ! Special treatment of magnetic field
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarInterleaveFormat1d(rvector1%NEQ, NVAR1D,&
                  'magneticfield_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat1d(rvector2%NEQ, NVAR1D,&
                  'magneticfield_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat1d(rvector3%NEQ, NVAR1D,&
                  'magneticfield_z', Ddata, p_Ddata3)

            case (NDIM2D)
              call mhd_getVarInterleaveFormat2d(rvector1%NEQ, NVAR2D,&
                  'magneticfield_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat2d(rvector2%NEQ, NVAR2D,&
                  'magneticfield_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat2d(rvector3%NEQ, NVAR2D,&
                  'magneticfield_z', Ddata, p_Ddata3)

            case (NDIM3D)
              call mhd_getVarInterleaveFormat3d(rvector1%NEQ, NVAR3D,&
                  'magneticfield_x', Ddata, p_Ddata1)
              call mhd_getVarInterleaveFormat3d(rvector2%NEQ, NVAR3D,&
                  'magneticfield_y', Ddata, p_Ddata2)
              call mhd_getVarInterleaveFormat3d(rvector3%NEQ, NVAR3D,&
                  'magneticfield_z', Ddata, p_Ddata3)
            end select

            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, 'magneticfield',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata2, 'magneticfield',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata3, 'magneticfield',&
                  density_ref, velocity_ref, length_ref)
            end if
            
            call ucd_addVarVertBasedVec(rexport, 'magneticfield'//csuffix,&
                p_Ddata1, p_Ddata2, p_Ddata3)

            if (btracers) then
              call ucd_addTracerVariable(rexport, 'magneticfield_x'//csuffix,&
                  p_Ddata1)
              call ucd_addTracerVariable(rexport, 'magneticfield_y'//csuffix,&
                  p_Ddata2)
              call ucd_addTracerVariable(rexport, 'magneticfield_z'//csuffix,&
                  p_Ddata3)
            end if
                       
          else

            ! Standard treatment for scalar quantity
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarInterleaveFormat1d(rvector1%NEQ,  NVAR1D,&
                  sucdvariable, Ddata, p_Ddata1)
            case (NDIM2D)
              call mhd_getVarInterleaveFormat2d(rvector1%NEQ,  NVAR2D,&
                  sucdvariable, Ddata, p_Ddata1)
            case (NDIM3D)
              call mhd_getVarInterleaveFormat3d(rvector1%NEQ,  NVAR3D,&
                  sucdvariable, Ddata, p_Ddata1)
            end select
            
            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, sucdvariable,&
                  density_ref, velocity_ref, length_ref)
            end if
            
            call ucd_addVariableVertexBased(rexport, sucdvariable//csuffix,&
                UCD_VAR_STANDARD, p_Ddata1)

            if (btracers) then
              call ucd_addTracerVariable(rexport, sucdvariable//csuffix,&
                  p_Ddata1)
            end if
            
          end if
        end do
        
      case (SYSTEM_BLOCKFORMAT)

        ! Loop over all variables
        do ivariable = 1, nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sucdvariable', sucdvariable, isubstring=ivariable)
          
          if (trim(sucdvariable) .eq. 'velocity') then

            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  'velocity_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat1d(rvector2%NEQ, NVAR1D,&
                  'velocity_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat1d(rvector3%NEQ, NVAR1D,&
                  'velocity_z', Ddata, p_Ddata3)
              
            case (NDIM2D)
              call mhd_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  'velocity_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat2d(rvector2%NEQ, NVAR2D,&
                  'velocity_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat2d(rvector3%NEQ, NVAR2D,&
                  'velocity_z', Ddata, p_Ddata3)

            case (NDIM3D)
              call mhd_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  'velocity_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat3d(rvector2%NEQ, NVAR3D,&
                  'velocity_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat3d(rvector3%NEQ, NVAR3D,&
                  'velocity_z', Ddata, p_Ddata3)
            end select
            
            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, 'velocity',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata2, 'velocity',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata3, 'velocity',&
                  density_ref, velocity_ref, length_ref)
            end if
            
            call ucd_addVarVertBasedVec(rexport, 'velocity'//csuffix,&
                UCD_VAR_VELOCITY, p_Ddata1, p_Ddata2, p_Ddata3)
            
            if (btracers) then
              call ucd_addTracerVariable(rexport, 'velocity_x'//csuffix,&
                  p_Ddata1)
              call ucd_addTracerVariable(rexport, 'velocity_y'//csuffix,&
                  p_Ddata2)
              call ucd_addTracerVariable(rexport, 'velocity_z'//csuffix,&
                  p_Ddata3)
            end if
            
          elseif (trim(sucdvariable) .eq. 'momentum') then

            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  'momentum_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat1d(rvector2%NEQ, NVAR1D,&
                  'momentum_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat1d(rvector3%NEQ, NVAR1D,&
                  'momentum_z', Ddata, p_Ddata3)

            case (NDIM2D)
              call mhd_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  'momentum_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat2d(rvector2%NEQ, NVAR2D,&
                  'momentum_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat2d(rvector3%NEQ, NVAR2D,&
                  'momentum_z', Ddata, p_Ddata3)

            case (NDIM3D)
              call mhd_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  'momentum_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat3d(rvector2%NEQ, NVAR3D,&
                  'momentum_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat3d(rvector3%NEQ, NVAR3D,&
                  'momentum_z', Ddata, p_Ddata3)
            end select
            
            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, 'momentum',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata2, 'momentum',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata3, 'momentum',&
                  density_ref, velocity_ref, length_ref)
            end if
            
            call ucd_addVarVertBasedVec(rexport, 'momentum'//csuffix,&
                p_Ddata1, p_Ddata2, p_Ddata3)

            if (btracers) then
              call ucd_addTracerVariable(rexport, 'momentum_x'//csuffix,&
                  p_Ddata1)
              call ucd_addTracerVariable(rexport, 'momentum_y'//csuffix,&
                  p_Ddata2)
              call ucd_addTracerVariable(rexport, 'momentum_z'//csuffix,&
                  p_Ddata3)
            end if
            
          elseif (trim(sucdvariable) .eq. 'magneticfield') then
            
            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  'magneticfield_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat1d(rvector2%NEQ, NVAR1D,&
                  'magneticfield_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat1d(rvector3%NEQ, NVAR1D,&
                  'magneticfield_z', Ddata, p_Ddata3)

            case (NDIM2D)
              call mhd_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  'magneticfield_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat2d(rvector2%NEQ, NVAR2D,&
                  'magneticfield_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat2d(rvector3%NEQ, NVAR2D,&
                  'magneticfield_z', Ddata, p_Ddata3)

            case (NDIM3D)
              call mhd_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  'magneticfield_x', Ddata, p_Ddata1)
              call mhd_getVarBlockFormat3d(rvector2%NEQ, NVAR3D,&
                  'magneticfield_y', Ddata, p_Ddata2)
              call mhd_getVarBlockFormat3d(rvector3%NEQ, NVAR3D,&
                  'magneticfield_z', Ddata, p_Ddata3)
            end select
            
            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, 'magneticfield',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata2, 'magneticfield',&
                  density_ref, velocity_ref, length_ref)
              call mhd_convertVariableDim(p_Ddata3, 'magneticfield',&
                  density_ref, velocity_ref, length_ref)
            end if
            
            call ucd_addVarVertBasedVec(rexport, 'magneticfield'//csuffix,&
                p_Ddata1, p_Ddata2, p_Ddata3)

            if (btracers) then
              call ucd_addTracerVariable(rexport, 'magneticfield_x'//csuffix,&
                  p_Ddata1)
              call ucd_addTracerVariable(rexport, 'magneticfield_y'//csuffix,&
                  p_Ddata2)
              call ucd_addTracerVariable(rexport, 'magneticfield_z'//csuffix,&
                  p_Ddata3)
            end if
            
          else
            
            ! Standard treatment for scalar quantity
            select case(ndim)
            case (NDIM1D)
              call mhd_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  sucdvariable, Ddata, p_Ddata1)
            case (NDIM2D)
              call mhd_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  sucdvariable, Ddata, p_Ddata1)
            case (NDIM3D)
              call mhd_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  sucdvariable, Ddata, p_Ddata1)
            end select
            
            if (bconvert) then
              call mhd_convertVariableDim(p_Ddata1, sucdvariable,&
                  density_ref, velocity_ref, length_ref)
            end if

            call ucd_addVariableVertexBased(rexport, sucdvariable//csuffix,&
                UCD_VAR_STANDARD, p_Ddata1)
            
            if (btracers) then
              call ucd_addTracerVariable(rexport, sucdvariable//csuffix,&
                  p_Ddata1)
            end if

          end if
        end do

      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_outputSolution')
        call sys_halt()
      end select

      ! Release temporal memory
      call lsyssc_releaseVector(rvector1)
      call lsyssc_releaseVector(rvector2)
      call lsyssc_releaseVector(rvector3)
      
    end subroutine outputSolution

  end subroutine mhd_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine mhd_outputStatistics(rtimerTotal, ssectionName, rcollection)

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
    real(DP) :: dfraction


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

    dfraction = 100.0_DP/rtimerTotal%delapsedReal

    call output_line('                                Real time (seconds) CPU time (seconds)       Frac')
    call output_line('Time for computing solution   : '//&
                     trim(sys_sdEL(p_rtimerSolution%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerSolution%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerSolution%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerSolution%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerSolution%delapsedCPU/&
                     (p_rtimerSolution%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for mesh adaptivity      : '//&
                     trim(sys_sdEL(p_rtimerAdaptation%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAdaptation%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAdaptation%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAdaptation%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAdaptation%delapsedCPU/&
                     (p_rtimerAdaptation%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for error estimation     : '//&
                     trim(sys_sdEL(p_rtimerErrorEstimation%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerErrorEstimation%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerErrorEstimation%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerErrorEstimation%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerErrorEstimation%delapsedCPU/&
                     (p_rtimerErrorEstimation%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for triangulation        : '//&
                     trim(sys_sdEL(p_rtimerTriangulation%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerTriangulation%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerTriangulation%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerTriangulation%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerTriangulation%delapsedCPU/&
                     (p_rtimerTriangulation%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for coefficient assembly : '//&
                     trim(sys_sdEL(p_rtimerAssemblyCoeff%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyCoeff%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAssemblyCoeff%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyCoeff%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAssemblyCoeff%delapsedCPU/&
                     (p_rtimerAssemblyCoeff%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for matrix assembly      : '//&
                     trim(sys_sdEL(p_rtimerAssemblyMatrix%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyMatrix%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAssemblyMatrix%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyMatrix%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAssemblyMatrix%delapsedCPU/&
                     (p_rtimerAssemblyMatrix%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for vector assembly      : '//&
                     trim(sys_sdEL(p_rtimerAssemblyVector%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyVector%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAssemblyVector%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyVector%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAssemblyVector%delapsedCPU/&
                     (p_rtimerAssemblyVector%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for pre-/post-processing : '//&
                     trim(sys_sdEL(p_rtimerPrePostprocess%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerPrePostprocess%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerPrePostprocess%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerPrePostprocess%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerPrePostprocess%delapsedCPU/&
                     (p_rtimerPrePostprocess%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_lbrk()
    call output_line('Time for total simulation     : '//&
                     trim(sys_sdEL(rtimerTotal%delapsedReal,3))//'           '//&
                     trim(sys_sdEL(rtimerTotal%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*rtimerTotal%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(rtimerTotal%delapsedCPU/&
                     (rtimerTotal%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_lbrk()

  end subroutine mhd_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine mhd_projectSolution(rsourceVector, rdestVector)
    
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
      call lsyssc_lumpMatrix(rmatrix1, LSYSSC_LUMP_DIAG)
      call lsyssc_lumpMatrix(rmatrix2, LSYSSC_LUMP_DIAG)

      ! Set the number of the scalar subvector to the collection structure
      rcollection%IquickAccess(2) = iblock
      
      ! Assemble the linear form for the scalar subvector
      call linf_buildVectorScalar(rform, .true.,&
          rdestVector%RvectorBlock(iblock),&
          fcoeff_buildVectorSc_sim=mhd_coeffVectorFE,&
          rcollection=rcollection)

      ! Compute the lumped L2-projection
      call lsyssc_invertedDiagMatVec(rmatrix2, rdestVector%RvectorBlock(iblock),&
          1.0_DP, rdestVector%RvectorBlock(iblock))
      
      ! Compute density-mass
      call lsyssc_duplicateVector(rsourceVector%RvectorBlock(iblock), rvector,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
      call lsyssc_matVec(rmatrix1, rsourceVector%RvectorBlock(iblock),&
          rvector, 1.0_DP, 0.0_DP)
      call lsyssc_getbase_double(rvector, p_Ddata)
      dmass1 = sum(p_Ddata)
      call lsyssc_releaseVector(rvector)
      
      ! Compute density-mass
      call lsyssc_duplicateVector(rdestVector%RvectorBlock(iblock), rvector,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
      call lsyssc_matVec(rmatrix2, rdestVector%RvectorBlock(iblock),&
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

  end subroutine mhd_projectSolution

  !*****************************************************************************

!<subroutine>

  subroutine mhd_outputFeSolution(rparlist, ssectionName,&
      rproblemLevel, rsolution)

!<description>
    ! This subroutine exports the FE-solution vector in RAW format
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution
!</input>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: sfesolution
    character(LEN=SYS_STRLEN) :: sfevariable

    ! persistent variable
    integer, save :: ifilenumber = 1

    ! local variables
    type(t_vectorScalar) :: rvector
    real(DP), dimension(:), pointer :: p_Ddata,p_DdataVariable
    integer :: isize,isystemFormat,cf,ndim,ivariable,nvariable

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName,&
                                'output', soutputName)
    call parlst_getvalue_string(rparlist, trim(soutputName),&
                                'sfesolution', sfesolution,'')
    if (trim(adjustl(sfesolution)) .eq. '') return

    call parlst_getvalue_int(rparlist, ssectionName,&
                             'isystemformat', isystemformat)

    ! Get number of variables to be written
    nvariable = max(1,&
        parlst_querysubstrings(rparlist,&
        trim(soutputName), 'sfevariable'))
    if (nvariable .eq. 0) return

    ! Set pointer
    call lsysbl_getbase_double(rsolution, p_Ddata)

    ! Calculate data
    isize = size(p_Ddata)/mhd_getNVAR(rproblemLevel)
    ndim  = rproblemLevel%rtriangulation%ndim

    ! Create auxiliary memory
    call lsyssc_createVector(rvector, isize, .false.)
    call lsyssc_getbase_double(rvector, p_DdataVariable)

    ! Open file for writing
    call io_openFileForWriting(trim(adjustl(sfesolution))//&
        '.'//trim(sys_si0(ifilenumber,5))//'.raw', cf,&
        SYS_REPLACE, bformatted=.false.)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Write individual variables to file
    select case(isystemFormat)
    case(SYSTEM_INTERLEAVEFORMAT)

      select case(ndim)
      case (NDIM1D)
        ! Loop over all variables
        do ivariable =1,nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sfevariable', sfevariable, isubstring=ivariable)
          
          ! Get variable data
          call mhd_getVarInterleaveFormat1d(rvector%NEQ,  NVAR1D,&
              trim(adjustl(sfevariable)), p_Ddata, p_DdataVariable)
          write(cf) isize,trim(adjustl(sfevariable))
          call vecio_writeArray(p_DdataVariable, cf, '')
        end do

      case (NDIM2D)
        ! Loop over all variables
        do ivariable =1,nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sfevariable', sfevariable, isubstring=ivariable)
          
          ! Get variable data
          call mhd_getVarInterleaveFormat2d(rvector%NEQ,  NVAR2D,&
              trim(adjustl(sfevariable)), p_Ddata, p_DdataVariable)
          write(cf) isize,trim(adjustl(sfevariable))
          call vecio_writeArray(p_DdataVariable, cf, '')
        end do

      case (NDIM3D)
        ! Loop over all variables
        do ivariable =1,nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sfevariable', sfevariable, isubstring=ivariable)
          
          ! Get variable data
          call mhd_getVarInterleaveFormat3d(rvector%NEQ,  NVAR3D,&
              trim(adjustl(sfevariable)), p_Ddata, p_DdataVariable)
          write(cf) isize,trim(adjustl(sfevariable))
          call vecio_writeArray(p_DdataVariable, cf, '')
        end do
      end select
      
    case (SYSTEM_BLOCKFORMAT)

      select case(ndim)
      case (NDIM1D)
        ! Loop over all variables
        do ivariable =1,nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sfevariable', sfevariable, isubstring=ivariable)
          
          ! Get variable data
          call mhd_getVarBlockFormat1d(rvector%NEQ,  NVAR1D,&
              trim(adjustl(sfevariable)), p_Ddata, p_DdataVariable)
          write(cf) isize,trim(adjustl(sfevariable))
          call vecio_writeArray(p_DdataVariable, cf, '')
        end do

      case (NDIM2D)
        ! Loop over all variables
        do ivariable =1,nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sfevariable', sfevariable, isubstring=ivariable)
          
          ! Get variable data
          call mhd_getVarBlockFormat2d(rvector%NEQ,  NVAR2D,&
              trim(adjustl(sfevariable)), p_Ddata, p_DdataVariable)
          write(cf) isize,trim(adjustl(sfevariable))
          call vecio_writeArray(p_DdataVariable, cf, '')
        end do

      case (NDIM3D)
        ! Loop over all variables
        do ivariable =1,nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sfevariable', sfevariable, isubstring=ivariable)
          
          ! Get variable data
          call mhd_getVarBlockFormat3d(rvector%NEQ,  NVAR3D,&
              trim(adjustl(sfevariable)), p_Ddata, p_DdataVariable)
          write(cf) isize,trim(adjustl(sfevariable))
          call vecio_writeArray(p_DdataVariable, cf, '')
        end do
      end select

    case default
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_outputFeSolution')
      call sys_halt()
    end select
    
    ! Release temporal memory
    call lsyssc_releaseVector(rvector)

  end subroutine mhd_outputFeSolution

end module mhd_postprocessing
