!##############################################################################
!# ****************************************************************************
!# <name> optcanalysis </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to perform error analysis for the solution
!# of the optimal control problem. Here, routines can be found that calculate
!# the functional, which is to be minimised by the theory:
!#
!# $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alpha/2||u||^2 $$
!#
!# The following routines can be found here:
!#
!# 2.) optcana_nonstatFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
!#
!# 3.) optcana_analyticalError
!#     -> Calculates the error ||y-y0|| of a given function y0 to an
!#        analytically given reference function y0.
!# </purpose>
!##############################################################################

module optcanalysis

  use fsystem
  use storage
  use genoutput
  
  use mprimitives
  use derivatives
  use boundary
  use cubature
  use triangulation
  use spatialdiscretisation
  use timediscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use pprocerror
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use collection
  
  use spacetimevectors
  use analyticsolution
  
  use constantsdiscretisation
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use assemblytemplates
  
  use structuresoperatorasm
  use spacematvecassembly
  
  use kktsystemspaces
  use kktsystem
  use spacesolver
  
  use user_callback

  implicit none
  
  private
  
!<types>

!<typeblock>

  ! Intermediate data needed during the computation of values concerning
  ! the functional J(.).
  type t_nonstFuncData
  
    ! Pointer to the underlying KKT system data.
    type(t_kktSystem), pointer :: p_rkktSystem => null()
  
    ! H1/2 control after the application of the Poincare-Steklow operator.
    type(t_controlSpace) :: rcontrolH12BdcPC
  
  end type
  
!</typeblock>

  public :: t_nonstFuncData

!</types>
  
  !public :: optcana_stationaryFunctional
  public :: optcana_prepareEval
  public :: optcana_doneEval
  public :: optcana_nonstatFunctional
  public :: optcana_nonstatFunctionalAtTime
  !public :: optcana_analyticalError
  
contains

  !****************************************************************************

!<subroutine>

  subroutine optcana_prepareEval (rnonstFuncData,rkktSystem,rkktSubsolverSet)

!<description>  
  ! Initialises the evaluation of the functional.
  ! Computes intermediate data.
!</description>

!<input>
  ! Underlying solution of the KKT system.
  type(t_kktSystem), intent(inout), target :: rkktSystem

  ! Subsolver set used for solving subproblems in the KKT system.
  type(t_kktSubsolverSet), intent(inout) :: rkktSubsolverSet
!</input>

!<output>
  ! Intermediate data that encapsules all information for the
  ! computation of J(.).
  type(t_nonstFuncData), intent(out) :: rnonstFuncData
!</output>    

!<subroutine>

    ! local variables
    type(t_spaceslSolverStat) :: rstatLocal

    ! Save some references.
    rnonstFuncData%p_rkktSystem => rkktSystem

    ! Initialises an additional control vector
    call kktsp_initControlVector (rnonstFuncData%rcontrolH12BdcPC,&
        rkktsystem%p_rcontrol%p_rvector%p_rspaceDiscr,&
        rkktsystem%p_rcontrol%p_rvector%p_rtimeDiscr)

    ! Apply the PC-Steklow operator to all controls in time.
    call kkt_controlApplyPCSteklow (rkktsystem,rnonstFuncData%rcontrolH12BdcPC,&
        rkktSubsolverSet,rstatLocal)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine optcana_doneEval (rnonstFuncData)

!<description>  
  ! Cleans up the evaluation of the functional J(.).
!</description>

!<inputoutput>
  ! Intermediate data that encapsules all information for the
  ! computation of J(.).
  type(t_nonstFuncData), intent(inout) :: rnonstFuncData
!</inputoutput>

!<subroutine>

    ! Release the temporary control.
    call kktsp_doneControlVector(rnonstFuncData%rcontrolH12BdcPC)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine optcana_fcalc_diffToTarget(Dintvalues,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the term ||y-z||^2_L2 at a point in time.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), optional, target :: rcollection
!</input>

!<output>
    ! Returns the values of the integral
    real(DP), dimension(:), intent(out) :: Dintvalues
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dx,dy, dtime
    real(DP) :: dminx,dmaxx,dminy,dmaxy
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dy1,p_Dy2,p_Dz1,p_Dz2
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    type(p_t_spacetimeOpAsmAnalyticData), target :: rp_ranalyticData
    type(t_spacetimeOpAsmAnalyticData), pointer :: p_ranalyticData
    real(DP), dimension(:), pointer :: p_DobservationArea
    
    ! Calcel if no FEM function is given.
    if (revalVectors%ncount .eq. 0) return

    ! Skip interleaved vectors.
    if (revalVectors%p_RvectorData(1)%bisInterleaved) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Current time
    dtime = rcollection%DquickAccess(5)

    ! Get the analytic data from the integer array
    rp_ranalyticData = transfer(rcollection%IquickAccess(:),rp_ranalyticData)
    p_ranalyticData => rp_ranalyticData%p_rdata

    ! Domain of interest
    dminx = -SYS_MAXREAL_DP
    dminy = -SYS_MAXREAL_DP
    dmaxx = SYS_MAXREAL_DP
    dmaxy = SYS_MAXREAL_DP
    
    p_DobservationArea => &
        p_ranalyticData%p_rsettingsOptControl%p_DobservationArea
    
    if (associated(p_DobservationArea)) then
      ! Put the observation area into the collection
      dminx = p_DobservationArea(1)
      dminy = p_DobservationArea(2)
      dmaxx = p_DobservationArea(3)
      dmaxy = p_DobservationArea(4)
    end if

    ! Which equation do we have?    
    select case (p_ranalyticData%p_rphysics%cequation)

    ! *************************************************************
    ! Stokes/Navier Stokes.
    ! *************************************************************
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dy1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)
      p_Dy2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC2D)
      p_Dz1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_FUNC2D)
      p_Dz2 => revalVectors%p_RvectorData(4)%p_Ddata(:,:,DER_FUNC2D)
      
      ! Evaluate the target function
      call getTargetData (1,p_Dz1)
      call getTargetData (2,p_Dz2)

      ! Calculate the integral...
      Dintvalues(:) = 0.0_DP
      
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          ! Calculate the error only in the observation area
          if ( (dx .ge. dminx) .and. &
               (dy .ge. dminy) .and. &
               (dx .le. dmaxx) .and. &
               (dy .le. dmaxy) ) then
               
            ! Multiply the values by the cubature weight and sum up
            ! into the (squared) L2 error:
            Dintvalues(1) = Dintvalues(1) + p_DcubWeight(icubp,iel) * &
                ( ( p_Dy1(icubp,iel)-p_Dz1(icubp,iel) )**2 + &
                  ( p_Dy2(icubp,iel)-p_Dz2(icubp,iel) )**2 )
                  
          end if
            
        end do ! icubp
      
      end do ! iel

    ! *************************************************************
    ! Heat equation
    ! *************************************************************
    case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dy1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)
      p_Dz1 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC2D)
    
      ! Evaluate the target function
      call getTargetData (1,p_Dz1)

      ! Calculate the integral...
      Dintvalues(:) = 0.0_DP
      
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Calculate the error only in the observation area
          if ( (dx .ge. dminx) .and. &
               (dy .ge. dminy) .and. &
               (dx .le. dmaxx) .and. &
               (dy .le. dmaxy) ) then

            ! Multiply the values by the cubature weight and sum up
            ! into the (squared) L2 error:
            Dintvalues(1) = Dintvalues(1) + p_DcubWeight(icubp,iel) * &
                  ( p_Dy1(icubp,iel)-p_Dz1(icubp,iel) )**2
                  
          end if
            
        end do ! icubp
      
      end do ! iel

    end select

  contains
  
    subroutine getTargetData (icomp,Ddata)
    
    ! Calculate component icomp of the target function and write the
    ! values in the cubature points to Ddata
    
    integer, intent(in) :: icomp
    real(DP), dimension(:,:), intent(inout) :: Ddata
    integer :: ierror

      ! Evaluate the analytic solution
      call ansol_evaluate (rcollection,"TARGET",icomp,&
          Ddata,npointsPerElement,nelements,p_Dpoints,rassemblyData%p_IelementList,ierror)
          
      select case (ierror)
      case (-1)
        ! Oops, that is the user defined function.
        call user_initCollectForVecAssembly (p_ranalyticData%p_rglobalData,&
            p_ranalyticData%p_rrhsDual%iid,icomp,dtime,rcollection%p_rnextCollection)
        call user_fct_Target (&
            Ddata, nelements, npointsPerElement, p_Dpoints, rcollection%p_rnextCollection)
        call user_doneCollectForVecAssembly (&
            p_ranalyticData%p_rglobalData,rcollection%p_rnextCollection)
        
        ierror = 0
      
      case (1)
        call output_line ("Error evaluating the target function.", &
                          OU_CLASS_ERROR,OU_MODE_STD,"optcana_evalFunction")
        call sys_halt()
      end select

    end subroutine
        
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine optcana_diffToTarget(dintvalue,dtime,rprimalSol,ranalyticData,rcubatureInfo)

!<description>  
    ! Calculates the term ||y(t)-z(t)||^2_L2 at a point t in time.
!</description>

!<input>
  ! Point in time where to evaluate.
  real(DP), intent(in) :: dtime
  
  ! Analytic data of the space-time problem
  type(t_spacetimeOpAsmAnalyticData), intent(inout), target :: ranalyticData
  
  ! Cubature information structure
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
!</input>

!<inputoutput>
  ! Solution vector
  type(t_primalSpace), intent(inout) :: rprimalSol
!</inputoutput>

!<output>
  ! Integral value
  real(DP), intent(out) :: dintvalue
!</output>

!</subroutine>

    ! local variables
    type(t_fev2Vectors) :: revalVectors
    type(t_vectorBlock), pointer :: p_rvector
    type(t_collection), target :: rcollection,rusercollection
    integer :: i,iindex
    type(p_t_spacetimeOpAsmAnalyticData), target :: rp_ranalyticData
    real(DP) :: dtimerel
    
    dintvalue = 0.0_DP
    
    ! Rescale the time.
    call mprim_linearRescale(dtime,&
        rprimalSol%p_rvector%p_rtimeDiscr%dtimeInit,rprimalSol%p_rvector%p_rtimeDiscr%dtimeMax,&
        0.0_DP,1.0_DP,dtimerel)
    dtimerel = min(max(0.0_DP,dtimerel),1.0_DP)

    ! Reserve memory in the pool for the solution
    iindex = -1
    call sptivec_getFreeBufferFromPool (rprimalSol%p_rvectorAccess,iindex,p_rvector)
    
    ! Calculate y(t)
    call sptivec_getTimestepDataByTime (rprimalSol%p_rvector, dtimerel, p_rvector)

    ! Append y to the set of functions to be evaluated.
    do i=1,p_rvector%nblocks
      call fev2_addVectorToEvalList(revalVectors,p_rvector%RvectorBlock(i),0)
    end do
    
    ! Add dummy vectors for z
    do i=1,p_rvector%nblocks
      call fev2_addDummyVectorToEvalList(revalVectors)
    end do
    
    ! Prepare the target function
    call collct_init(rcollection)
    rcollection%p_rnextCollection => rusercollection
    call ansol_prepareEval (ranalyticData%p_rsettingsOptControl%rtargetFunction,&
        rcollection,"TARGET",dtime)
        
    ! Put the analytic data to the integer array of the collection
    rp_ranalyticData%p_rdata => ranalyticData
    rcollection%IquickAccess(:) = transfer(rp_ranalyticData,rcollection%IquickAccess(:))
    
    ! Calculate the integral
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,optcana_fcalc_diffToTarget,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)
    
    ! Cleanup
    call ansol_doneEvalCollection (rcollection,"TARGET")
    call collct_done(rcollection)

    call fev2_releaseVectorList(revalVectors)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine optcana_controlNormDistL2(dintvalue,dtime,rnonstFuncData,&
      icompstart,ncomponents,rcubatureInfo)

!<description>  
    ! Calculates the term ||u(t)||^2_L2 at a point t in time.
    ! u(.) is a distributed control variable.
!</description>

!<input>
  ! Point in time where to evaluate.
  real(DP), intent(in) :: dtime
  
  ! Number of the first component in the control defining u
  integer, intent(in) :: icompstart

  ! Number of components in u
  integer, intent(in) :: ncomponents
  
  ! Cubature information structure
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
!</input>

!<inputoutput>
  ! Intermediate data that encapsules all information for the
  ! computation of J(.).
  type(t_nonstFuncData), intent(inout) :: rnonstFuncData
!</inputoutput>

!<output>
  ! Integral value
  real(DP), intent(out) :: dintvalue
!</output>

!</subroutine>

    ! local variables
    type(t_fev2Vectors) :: revalVectors
    type(t_vectorBlock), pointer :: p_rvector
    integer :: i
    
    dintvalue = 0.0_DP
    
    ! Calculate u(t).
    call kkt_getControlAtTime (rnonstFuncData%p_rkktsystem, dtime, p_rvector)

    ! Append u to the set of functions to be evaluated.
    do i=icompstart,icompstart+ncomponents-1
      call fev2_addVectorToEvalList(revalVectors,p_rvector%RvectorBlock(i),0)
    end do
    
    ! Calculate the integral
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,bma_fcalc_L2norm,&
        revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    
    ! Cleanup
    call fev2_releaseVectorList(revalVectors)
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_getL2BDControlValues (Dvalues,ibct,DpointPar,rbdControlVec,&
      rtriangulation,rboundary)
  
!<description>
  ! Calculates the values of a FEM function in a couple of points on
  ! the boundary, given by a boundary component and a set of parameter values.
  !
  ! WORKS ONLY FOR P2 ELEMENTS ON THE BOUNDARY !!!
!</description>

!<input>
  ! Boundary component
  integer, intent(in) :: ibct
  
  ! Parameter values of the points, in length parametrisation
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(in) :: DpointPar
  
  ! Underlying 2D Triangulation.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Underlying boundary definition
  type(t_boundary), intent(in) :: rboundary
  
  ! Control vector on the boundary
  type(t_vectorScalar), intent(in) :: rbdControlVec
!</input>

!<output>
  ! Values of the boudary control vector in the points
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    integer :: npointsPerElement,nelements,iel,ipt,iindex,ivt1,ivt2,nvtbd
    real(DP) :: d1,d2,d3,dpar,dstartpar,dendpar,dparrel
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    real(DP), dimension(:), pointer :: p_Ddata,p_DvertexParameterValue
    
    ! Getch some data.
    call storage_getbase_int2d(&
        rbdControlVec%p_rspatialDiscr%p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(&
        rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
    call lsyssc_getbase_double (rbdControlVec,p_Ddata)
    nvtbd = rbdControlVec%p_rspatialDiscr%p_rtriangulation%nvt
    
    npointsPerElement = ubound(Dvalues,1)
    nelements = ubound(Dvalues,2)

    ! Loop over all points    
    do iel = 1,nelements
      do ipt = 1,npointsPerElement
      
        ! Get the index of the edge in the 2D parametrisation
        call tria_searchBoundaryEdgePar2D(ibct, DpointPar(ipt,iel), &
            rtriangulation, rboundary, iindex, BDR_PAR_LENGTH)
            
        ! The edge index in 2D is the element index in the
        ! 1D parametrisation. Get the associated points.
        ivt1 = p_IverticesAtElement(1,iindex)
        ivt2 = p_IverticesAtElement(2,iindex)
        
        ! Get the start and end parameter value of the current edge.
        dpar = DpointPar(ipt,iel)
        dstartpar = p_DvertexParameterValue(iindex)
        if (iindex .lt. p_IboundaryCpIdx(ibct+1)-1) then
          dendpar = p_DvertexParameterValue(iindex+1)
        else
          dendpar = boundary_dgetMaxParVal(rboundary,ibct)
        end if
        
        ! Calculate the "relative" parameter value. Scale dpar to [-1,1]
        call mprim_linearRescale(dpar,dstartpar,dendpar,-1.0_DP,1.0_DP,dparrel)
        
        ! Get the values of the P2 element on the boundary.
        d1 = p_Ddata(ivt1)
        d2 = p_Ddata(iindex+nvtbd)
        d3 = p_Ddata(ivt2)
        
        ! Quadratic interpolation, calculate the value at dpar.
        call mprim_quadraticInterpolation (dparrel,d1,d2,d3,Dvalues(ipt,iel))
      
      end do
    end do  

  end subroutine

!******************************************************************************

  !<subroutine>

    subroutine ffunction_dirichletL2BdC (cderivative, rdiscretisation, &
                                   DpointsRef, Dpoints, ibct, DpointPar,&
                                   Ielements, Dvalues, rcollection)

    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection

  !<description>
    ! Routine to calculate the values of a 1D function on the boundary.
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

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

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
  
      ! local variables
      type(t_vectorBlock), pointer :: p_rvector
      integer :: icomponent
      
      p_rvector => rcollection%p_rvectorQuickAccess1
      icomponent = rcollection%IquickAccess(1)
      
      ! Calculate the values of the boundary function in the points.
      call optcana_getL2BDControlValues (Dvalues,ibct,DpointPar,&
          p_rvector%RvectorBlock(icomponent),&
          rdiscretisation%p_rtriangulation,rdiscretisation%p_rboundary)

    end subroutine

  !****************************************************************************

!<subroutine>

  subroutine optcana_controlNormBdCL2(dintvalue,dtime,rnonstFuncData,&
      icompstart,ncomponents)

!<description>  
    ! Calculates the term ||u(t)||^2_L2 at a point t in time.
    ! u(.) is a boundary control variable.
!</description>

!<input>
  ! Point in time where to evaluate.
  real(DP), intent(in) :: dtime
  
  ! Number of the first component in the control defining u
  integer, intent(in) :: icompstart

  ! Number of components in u
  integer, intent(in) :: ncomponents
!</input>

!<inputoutput>
  ! Intermediate data that encapsules all information for the
  ! computation of J(.).
  type(t_nonstFuncData), intent(inout) :: rnonstFuncData
!</inputoutput>

!<output>
  ! Integral value
  real(DP), intent(out) :: dintvalue
!</output>

!</subroutine>

    ! local variables
    type(t_fev2Vectors) :: revalVectors
    type(t_vectorBlock), pointer :: p_rvector
    integer :: i
    type(t_collection) :: rcollection
    real(DP) :: derror
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    
    p_rdiscretisation => &
        rnonstFuncData%p_rkktsystem%p_rprimalSol%p_rvector%p_rspaceDiscr%RspatialDiscr(1)
    
    dintvalue = 0.0_DP
    
    ! Calculate u(t).
    call kkt_getControlAtTime (rnonstFuncData%p_rkktsystem, dtime, p_rvector)

    ! Loop over the vector components to calculate
    do i=icompstart,icompstart+ncomponents-1
    
      ! Prepare the collection
      rcollection%p_rvectorQuickAccess1 => p_rvector
      rcollection%IquickAccess(1) = i
      
      ! Calculate the integral
      call pperr_scalarBoundary2D (PPERR_L2ERROR, CUB_G4_1D, derror,&
          ffunctionReference=ffunction_dirichletL2BdC, rcollection=rcollection, &
          rdiscretisation=p_rdiscretisation)
          
      dintvalue = dintvalue + derror**2
    
    end do

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_getH12BDControlValues (Dvalues,ibct,DpointPar,&
      rbdControlVec1,rbdControlVec2,rtriangulation,rboundary)
  
!<description>
  ! Calculates the values of a the function <u,Su> in a couple of points on
  ! the boundary, given by a boundary component and a set of parameter values.
  ! "Su" is the function which is returned by the Steklow-Poincare operator
  ! and muzst be specified in rbdControlVec2.
  !
  ! WORKS ONLY FOR P2 ELEMENTS ON THE BOUNDARY !!!
!</description>

!<input>
  ! Boundary component
  integer, intent(in) :: ibct
  
  ! Parameter values of the points, in length parametrisation
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(in) :: DpointPar
  
  ! Underlying 2D Triangulation.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Underlying boundary definition
  type(t_boundary), intent(in) :: rboundary
  
  ! Control vector on the boundary
  type(t_vectorScalar), intent(in) :: rbdControlVec1

  ! Control vector on the boundary after applying the Steklow-Poincare operator
  type(t_vectorScalar), intent(in) :: rbdControlVec2
!</input>

!<output>
  ! Values of the boudary control vector in the points
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
  
!</subroutine>

    ! local variables
    integer :: npointsPerElement,nelements,iel,ipt,iindex,ivt1,ivt2,nvtbd
    real(DP) :: d1,d2,d3,dpar,dstartpar,dendpar,dparrel,dvalue1,dvalue2
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_DvertexParameterValue
    
    ! Getch some data.
    call storage_getbase_int2d(&
        rbdControlVec1%p_rspatialDiscr%p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(&
        rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
    call lsyssc_getbase_double (rbdControlVec1,p_Ddata1)
    call lsyssc_getbase_double (rbdControlVec2,p_Ddata2)
    nvtbd = rbdControlVec1%p_rspatialDiscr%p_rtriangulation%nvt
    
    npointsPerElement = ubound(Dvalues,1)
    nelements = ubound(Dvalues,2)

    ! Loop over all points    
    do iel = 1,nelements
      do ipt = 1,npointsPerElement
      
        ! Get the index of the edge in the 2D parametrisation
        call tria_searchBoundaryEdgePar2D(ibct, DpointPar(ipt,iel), &
            rtriangulation, rboundary, iindex, BDR_PAR_LENGTH)
            
        ! The edge index in 2D is the element index in the
        ! 1D parametrisation. Get the associated points.
        ivt1 = p_IverticesAtElement(1,iindex)
        ivt2 = p_IverticesAtElement(2,iindex)
        
        ! Get the start and end parameter value of the current edge.
        dpar = DpointPar(ipt,iel)
        dstartpar = p_DvertexParameterValue(iindex)
        if (iindex .lt. p_IboundaryCpIdx(ibct+1)-1) then
          dendpar = p_DvertexParameterValue(iindex+1)
        else
          dendpar = boundary_dgetMaxParVal(rboundary,ibct)
        end if
        
        ! Calculate the "relative" parameter value. Scale dpar to [-1,1]
        call mprim_linearRescale(dpar,dstartpar,dendpar,-1.0_DP,1.0_DP,dparrel)
        
        ! Get the values of the P2 element on the boundary. Control.
        d1 = p_Ddata1(ivt1)
        d2 = p_Ddata1(iindex+nvtbd)
        d3 = p_Ddata1(ivt2)
        
        ! Quadratic interpolation, calculate the value at dpar.
        call mprim_quadraticInterpolation (dparrel,d1,d2,d3,dvalue1)

        ! Get the values of the P2 element on the boundary. Preconditioned control "Su".
        d1 = p_Ddata2(ivt1)
        d2 = p_Ddata2(iindex+nvtbd)
        d3 = p_Ddata2(ivt2)
        
        ! Quadratic interpolation, calculate the value at dpar.
        call mprim_quadraticInterpolation (dparrel,d1,d2,d3,dvalue2)
      
        ! Return <u,Su>
        Dvalues(ipt,iel) = dvalue1 * dvalue2
      end do
    end do  

  end subroutine

!******************************************************************************

  !<subroutine>

    subroutine ffunction_dirichletH12BdC (cderivative, rdiscretisation, &
                                   DpointsRef, Dpoints, ibct, DpointPar,&
                                   Ielements, Dvalues, rcollection)

    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection

  !<description>
    ! Routine to calculate the values of a 1D function on the boundary.
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

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

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
  
      ! local variables
      type(t_vectorBlock), pointer :: p_rvector1,p_rvector2
      integer :: icomponent
      
      p_rvector1 => rcollection%p_rvectorQuickAccess1
      p_rvector2 => rcollection%p_rvectorQuickAccess2
      icomponent = rcollection%IquickAccess(1)
      
      ! Calculate the values of the boundary function in the points.
      call optcana_getH12BDControlValues (Dvalues,ibct,DpointPar,&
          p_rvector1%RvectorBlock(icomponent),&
          p_rvector2%RvectorBlock(icomponent),&
          rdiscretisation%p_rtriangulation,rdiscretisation%p_rboundary)

    end subroutine

  !****************************************************************************

!<subroutine>

  subroutine optcana_controlNormBdCH12(dintvalue,dtime,rnonstFuncData,&
      icompstart,ncomponents)

!<description>  
    ! Calculates the term ||u(t)||^2_H12 at a point t in time.
    ! u(.) is a boundary control variable.
!</description>

!<input>
  ! Point in time where to evaluate.
  real(DP), intent(in) :: dtime
  
  ! Number of the first component in the control defining u
  integer, intent(in) :: icompstart

  ! Number of components in u
  integer, intent(in) :: ncomponents
!</input>

!<inputoutput>
  ! Intermediate data that encapsules all information for the
  ! computation of J(.).
  type(t_nonstFuncData), intent(inout) :: rnonstFuncData
!</inputoutput>

!<output>
  ! Integral value
  real(DP), intent(out) :: dintvalue
!</output>

!</subroutine>

    ! local variables
    type(t_fev2Vectors) :: revalVectors
    type(t_vectorBlock), pointer :: p_rvector1,p_rvector2
    integer :: i
    type(t_collection) :: rcollection
    real(DP) :: derror
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    
    p_rdiscretisation => &
        rnonstFuncData%p_rkktsystem%p_rprimalSol%p_rvector%p_rspaceDiscr%RspatialDiscr(1)
    
    dintvalue = 0.0_DP
    
    ! Calculate u(t).
    call kkt_getControlAtTime (rnonstFuncData%p_rkktsystem, dtime, p_rvector1)
    call kkt_getControlAtTime (rnonstFuncData%p_rkktsystem, dtime, p_rvector2, &
        rnonstFuncData%rcontrolH12BdcPC)

    ! Loop over the vector components to calculate
    do i=icompstart,icompstart+ncomponents-1
    
      ! Prepare the collection
      rcollection%p_rvectorQuickAccess1 => p_rvector1
      rcollection%p_rvectorQuickAccess2 => p_rvector2
      rcollection%IquickAccess(1) = i
      
      ! Calculate the integral.
      ! PPERR_MEANERROR just calculates "int(f)", and here is f=<u,Su>.
      call pperr_scalarBoundary2D (PPERR_MEANERROR, CUB_G4_1D, derror,&
          ffunctionReference=ffunction_dirichletH12BdC, rcollection=rcollection, &
          rdiscretisation=p_rdiscretisation)
          
      dintvalue = dintvalue + derror
    
    end do

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_nonstatFunctional (Derror, &
      ispacelevel, itimelevel, roperatorAsmHier, rnonstFuncData)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions
  ! user_ffunction_TargetX and user_ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
!</description>
  
!<input>
  ! Space-level corresponding to the solutions
  integer, intent(in) :: ispacelevel

  ! Time-level corresponding to the solutions
  integer, intent(in) :: itimelevel
  
  ! Hierarchy of space-time operators.
  type(t_spacetimeOpAsmHierarchy), intent(in) :: roperatorAsmHier
!</input>

!<inputoutput>
  ! Intermediate data that encapsules all information for the
  ! computation of J(.).
  type(t_nonstFuncData), intent(inout) :: rnonstFuncData
!</inputoutput>

!<output>
  ! Returns information about the error.
  ! Derror(1) = J(y,u).
  ! Derror(2) = ||y-z||_{L^2}.
  ! Derror(3) = ||y(T)-z(T)||_{L^2}.
  ! Derror(4) = ||u||_{L^2}.
  ! Derror(5) = ||u||_{L^2(Gamma_C)}.
  ! Derror(6) = (Su,u)                   (with S the Steklov-Poincare operator)
  real(DP), dimension(:), intent(out) :: Derror
!</output>
  
!</subroutine>
    
    ! local variables
    type(t_spacetimeOperatorAsm) :: roperatorAsm
    
    real(DP) :: dtime,dtimestart,dtimeend,dtstep,dval,dtheta
    
    integer :: istep,idoftime,icomp,ncomp
    real(DP),dimension(2) :: Derr
    type(t_collection), target :: rcollection,rlocalcoll
    type(t_vectorBlock), target :: rtempVector, rzeroVector
    real(dp), dimension(:), pointer :: p_DobservationArea
    real(dp), dimension(:), pointer :: p_Dx
    type(t_controlSpace) :: rcontrolPC
    type(t_dualSpace) :: rdualPC
    type(t_vectorBlock), pointer :: p_rcontrolVec,p_rcontrolPCVec
    type(t_spaceslSolverStat) :: rstatLocal
    integer :: ierror
    
    !type(t_sptiDirichletBCCBoundary) :: rdirichletBCC
    !type(t_bdRegionEntry), pointer :: p_rbdRegionIterator
    
    ! Get the corresponding operator assembly structure
    call stoh_getOpAsm_slvtlv (&
        roperatorAsm,roperatorAsmHier,ispacelevel,itimelevel)
        
    ! Initialise the collection for the assembly process with callback routines.
    ! This stores the simulation time in the collection and sets the
    ! current subvector z for the callback routines.
    call collct_init(rcollection)
    
    ! Clear the output
    Derror(1:6) = 0.0_DP

!    ! Assemble the dirichlet control boundary conditions
!    call stdbcc_createDirichletBCCBd (rsolution%p_rspaceDiscr,rsolution%p_rtimeDiscr,&
!        rdirichletBCC)
!    call stdbcc_assembleDirichletBCCBd (roptcBDC,rdirichletBCC,rglobalData)

    ! Timestepping technique?
    select case (roperatorAsm%p_rtimeDiscrPrimal%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = roperatorAsm%p_rtimeDiscrPrimal%dtheta

      ! Characteristics of the timestepping is taken from the 1st timestep
      call tdiscr_getTimestep(roperatorAsm%p_rtimeDiscrPrimal,1,dtstep=dtstep)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (roperatorAsm%p_rtimeDiscrPrimal%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getRhs_Primal")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! ---------------------------------------------------------
        ! Calculate ||y-z||^2
        ! ---------------------------------------------------------
        
        do idoftime = 1,rnonstFuncData%p_rkktsystem%p_rprimalSol%p_rvector%NEQtime
        
          ! Characteristics of the current timestep.
          call tdiscr_getTimestep(roperatorAsm%p_rtimeDiscrPrimal,idofTime-1,&
              dtimeend,dtstep,dtimestart)
              
          ! Calculate ||y-z||^2 in the endpoint of the time interval.
          call optcana_diffToTarget(dval,dtimeend,&
              rnonstFuncData%p_rkktsystem%p_rprimalSol,&
              roperatorAsm%p_ranalyticData,roperatorAsm%p_rasmTemplates%rcubatureInfoRHS)
              
          ! Sum up according to the summed trapezoidal rule
          if ((idoftime .eq. 1) .or. &
              (idoftime .eq. rnonstFuncData%p_rkktsystem%p_rprimalSol%p_rvector%NEQtime)) then
            Derror(2) = Derror(2) + 0.5_DP * dval * dtstep
          else
            Derror(2) = Derror(2) + dval * dtstep
          end if
        
          ! ---------------------------------------------------------
          ! Calculate ||y(T)-z(T)||^2
          ! ---------------------------------------------------------
          if ((idoftime .eq. &
               rnonstFuncData%p_rkktsystem%p_rprimalSol%p_rvector%NEQtime) .and. &
              (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%ddeltaC .gt. 0.0_DP)) then
            Derror(3) = dval
          end if
          
        end do

        ! ---------------------------------------------------------
        ! Calculate ||u||^2
        ! ---------------------------------------------------------
        icomp = 1

        ! ---------------------------------------------------------
        ! Distributed control
        ! ---------------------------------------------------------
        if (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

          ! Type of equation?
          select case (roperatorAsmHier%ranalyticData%p_rphysics%cequation)

          ! *************************************************************
          ! Stokes/Navier Stokes.
          ! *************************************************************
          case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
            ncomp = 2

          ! *************************************************************
          ! Heat equation
          ! *************************************************************
          case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
            ncomp = 1

          end select

          do idoftime = 1,rnonstFuncData%p_rkktsystem%p_rcontrol%p_rvector%NEQtime
          
            ! Characteristics of the current timestep.
            call tdiscr_getTimestep(roperatorAsm%p_rtimeDiscrPrimal,idofTime-1,&
                dtimeend,dtstep,dtimestart)
                
            ! The control lives at the point dtime in time.
            dtime = (1.0_DP-dtheta)*dtimestart + dtheta*dtimeend

            ! Calculate ||u||^2 at the point in time where u lives.
            ! Note that we pass dtimeend here as time. This is correct!
            ! dtimeend on the time scale of rcontrol corresponds to the time dtime
            ! as rcontrol is shifted by a half timestep!
            call optcana_controlNormDistL2(dval,dtimeend,&
                rnonstFuncData,icomp,ncomp,&
                roperatorAsm%p_rasmTemplates%rcubatureInfoRHS)

            ! Sum up according to the rectangular rule
            Derror(4) = Derror(4) + dval * dtstep
          
          end do
          
          icomp = icomp + ncomp
          
        end if

        ! ---------------------------------------------------------
        ! L2 boundary control
        ! ---------------------------------------------------------
        if (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaL2BdC .ge. 0.0_DP) then

          ! Type of equation?
          select case (roperatorAsmHier%ranalyticData%p_rphysics%cequation)

          ! *************************************************************
          ! Stokes/Navier Stokes.
          ! *************************************************************
          case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
            ncomp = 2

          ! *************************************************************
          ! Heat equation
          ! *************************************************************
          case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
            ncomp = 1

          end select

          do idoftime = 1,rnonstFuncData%p_rkktsystem%p_rcontrol%p_rvector%NEQtime
          
            ! Characteristics of the current timestep.
            call tdiscr_getTimestep(roperatorAsm%p_rtimeDiscrPrimal,idofTime-1,&
                dtimeend,dtstep,dtimestart)
                
            ! The control lives at the point dtime in time.
            dtime = (1.0_DP-dtheta)*dtimestart + dtheta*dtimeend

            ! Calculate ||u||^2 at the point in time where u lives.
            ! Note that we pass dtimeend here as time. This is correct!
            ! dtimeend on the time scale of rcontrol corresponds to the time dtime
            ! as rcontrol is shifted by a half timestep!
            call optcana_controlNormBdCL2(dval,dtimeend,rnonstFuncData,icomp,ncomp)

            ! Sum up according to the rectangular rule
            Derror(5) = Derror(5) + dval * dtstep
          
          end do
          
          icomp = icomp + ncomp
          
        end if
      
        ! ---------------------------------------------------------
        ! H1/2 boudary control
        ! ---------------------------------------------------------
        if (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaH12BdC .ge. 0.0_DP) then

          ! Type of equation?
          select case (roperatorAsmHier%ranalyticData%p_rphysics%cequation)

          ! *************************************************************
          ! Stokes/Navier Stokes.
          ! *************************************************************
          case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
            ncomp = 2

          ! *************************************************************
          ! Heat equation
          ! *************************************************************
          case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
            ncomp = 1

          end select
          
          ! Calculate the actual time integral.
          do idoftime = 1,rnonstFuncData%p_rkktsystem%p_rcontrol%p_rvector%NEQtime
          
            ! Characteristics of the current timestep.
            call tdiscr_getTimestep(roperatorAsm%p_rtimeDiscrPrimal,idofTime-1,&
                dtimeend,dtstep,dtimestart)
                
            ! The control lives at the point dtime in time.
            dtime = (1.0_DP-dtheta)*dtimestart + dtheta*dtimeend
            ! Actually not used... Below, dtimeend is used.

            ! Calculate ||u||^2 at the point in time where u lives.
            ! Note that we pass dtimeend here as time. This is correct!
            ! dtimeend on the time scale of rcontrol corresponds to the time dtime
            ! as rcontrol is shifted by a half timestep!
            call optcana_controlNormBdCH12(dval,dtimeend,rnonstFuncData,icomp,ncomp)

            ! Sum up according to the rectangular rule
            Derror(6) = Derror(6) + dval * dtstep
          
          end do
          
          icomp = icomp + ncomp
          
        end if

      case default      
        call output_line("Unknown timestep sub-scheme",&
            OU_CLASS_ERROR,OU_MODE_STD,"optcana_nonstatFunctional")
        call sys_halt()
      end select ! Tag

    case default      
      call output_line("Unknown timestep scheme",&
          OU_CLASS_ERROR,OU_MODE_STD,"optcana_nonstatFunctional")
      call sys_halt()
    end select ! Timestep scheme

    ! Clean up the collection
    call collct_done(rcollection)
      
    ! Calculate J(.)
    Derror(1) = 0.5_DP*Derror(2)  &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaEndTimeC * Derror(3)  &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaDistC * Derror(4) &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaL2BdC * Derror(5) &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaH12BdC * Derror(6)
              
    ! Take some square roots to calculate the actual values.
    Derror(2:) = sqrt(Derror(2:))
    

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_nonstatFunctionalAtTime (Derror, dtime,&
      ispacelevel, itimelevel, roperatorAsmHier, rnonstFuncData)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions
  ! user_ffunction_TargetX and user_ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
  !
  ! The subroutine returns the values of the space integrals at a special
  ! point in time. THere is no time integratino applied.
!</description>
  
!<input>
  ! Point in time where to evaluate the values of the functional.
  real(DP), intent(in) :: dtime

  ! Space-level corresponding to the solutions
  integer, intent(in) :: ispacelevel

  ! Time-level corresponding to the solutions
  integer, intent(in) :: itimelevel
  
  ! Hierarchy of space-time operators.
  type(t_spacetimeOpAsmHierarchy), intent(in) :: roperatorAsmHier
!</input>

!<inputoutput>
  ! Intermediate data that encapsules all information for the
  ! computation of J(.).
  type(t_nonstFuncData), intent(inout) :: rnonstFuncData
!</inputoutput>

!<output>
  ! Returns information about the error.
  ! Derror(1) = J(y,u).
  ! Derror(2) = ||y-z||_{L^2}.
  ! Derror(3) = ||y(T)-z(T)||_{L^2}.
  ! Derror(4) = ||u||_{L^2}.
  ! Derror(5) = ||u||_{L^2(Gamma_C)}.
  ! Derror(6) = ||u||_{H^1/2(Gamma_C)}.
  real(DP), dimension(:), intent(out) :: Derror
!</output>
  
!</subroutine>
    
    ! local variables
    type(t_spacetimeOperatorAsm) :: roperatorAsm
    
    real(DP) :: dval
    
    integer :: istep,icomp,ncomp
    real(DP),dimension(2) :: Derr
    type(t_collection), target :: rcollection,rlocalcoll
    type(t_vectorBlock), target :: rtempVector, rzeroVector
    real(dp), dimension(:), pointer :: p_DobservationArea
    real(dp), dimension(:), pointer :: p_Dx
    !type(t_sptiDirichletBCCBoundary) :: rdirichletBCC
    !type(t_bdRegionEntry), pointer :: p_rbdRegionIterator
    
    ! Get the corresponding operator assembly structure
    call stoh_getOpAsm_slvtlv (&
        roperatorAsm,roperatorAsmHier,ispacelevel,itimelevel)
        
    ! Initialise the collection for the assembly process with callback routines.
    ! This stores the simulation time in the collection and sets the
    ! current subvector z for the callback routines.
    call collct_init(rcollection)
    
    ! Clear the output
    Derror(1:6) = 0.0_DP

    ! ---------------------------------------------------------
    ! Calculate ||y(t)-z(t)||^2
    ! ---------------------------------------------------------
    
    ! Calculate ||y-z||^2 in the endpoint of the time interval.
    call optcana_diffToTarget(dval,dtime,rnonstFuncData%p_rkktsystem%p_rprimalSol,&
        roperatorAsm%p_ranalyticData,roperatorAsm%p_rasmTemplates%rcubatureInfoRHS)
    
    Derror(2) = dval
    
    ! ---------------------------------------------------------
    ! Calculate ||y(T)-z(T)||^2
    ! ---------------------------------------------------------
    if ((dtime .eq. rnonstFuncData%p_rkktsystem%p_rprimalSol%p_rvector%p_rtimeDiscr%dtimeMax) .and. &
        (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%ddeltaC .gt. 0.0_DP)) then
      Derror(3) = dval
    end if
      
    ! ---------------------------------------------------------
    ! Calculate ||u||^2
    ! ---------------------------------------------------------
    icomp = 1

    ! ---------------------------------------------------------
    ! Distributed control
    ! ---------------------------------------------------------
    if (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

      ! Type of equation?
      select case (roperatorAsmHier%ranalyticData%p_rphysics%cequation)

      ! *************************************************************
      ! Stokes/Navier Stokes.
      ! *************************************************************
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
        ncomp = 2

      ! *************************************************************
      ! Heat equation
      ! *************************************************************
      case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
        ncomp = 1

      end select

      ! Calculate ||u||^2 at the point in time where u lives.
      ! Note that we pass dtimeend here as time. This is correct!
      ! dtimeend on the time scale of rcontrol corresponds to the time dtime
      ! as rcontrol is shifted by a half timestep!
      call optcana_controlNormDistL2(dval,dtime,&
          rnonstFuncData,icomp,ncomp,&
          roperatorAsm%p_rasmTemplates%rcubatureInfoRHS)

      Derror(4) = dval
      
      icomp = icomp + ncomp
      
    end if

    ! ---------------------------------------------------------
    ! L2 boundary control
    ! ---------------------------------------------------------
    if (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaL2BdC .ge. 0.0_DP) then

      ! Type of equation?
      select case (roperatorAsmHier%ranalyticData%p_rphysics%cequation)

      ! *************************************************************
      ! Stokes/Navier Stokes.
      ! *************************************************************
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
        ncomp = 2

      ! *************************************************************
      ! Heat equation
      ! *************************************************************
      case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
        ncomp = 1

      end select

      ! Calculate ||u||^2 at the point in time where u lives.
      ! Note that we pass dtimeend here as time. This is correct!
      ! dtimeend on the time scale of rcontrol corresponds to the time dtime
      ! as rcontrol is shifted by a half timestep!
      call optcana_controlNormBdCL2(dval,dtime,rnonstFuncData,icomp,ncomp)

      Derror(5) = dval
      
      icomp = icomp + ncomp
      
    end if

    ! ---------------------------------------------------------
    ! L2 boundary control
    ! ---------------------------------------------------------
    if (roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaH12BdC .ge. 0.0_DP) then

      ! Type of equation?
      select case (roperatorAsmHier%ranalyticData%p_rphysics%cequation)

      ! *************************************************************
      ! Stokes/Navier Stokes.
      ! *************************************************************
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
        ncomp = 2

      ! *************************************************************
      ! Heat equation
      ! *************************************************************
      case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
        ncomp = 1

      end select

      ! Calculate ||u||^2 at the point in time where u lives.
      ! Note that we pass dtimeend here as time. This is correct!
      ! dtimeend on the time scale of rcontrol corresponds to the time dtime
      ! as rcontrol is shifted by a half timestep!
      call optcana_controlNormBdCH12(dval,dtime,rnonstFuncData,icomp,ncomp)

      Derror(6) = dval
      
      icomp = icomp + ncomp
      
    end if

    ! Clean up the collection
    call collct_done(rcollection)
      
    ! Calculate J(.)
    Derror(1) = 0.5_DP*Derror(2)  &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaEndTimeC * Derror(3)  &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaDistC * Derror(4) &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaL2BdC * Derror(5) &
              + 0.5_DP*roperatorAsmHier%ranalyticData%p_rsettingsOptControl%dalphaL2BdC * Derror(6)
              
    ! Take some square roots to calculate the actual values.
    Derror(2:) = sqrt(Derror(2:))

  end subroutine

!!******************************************************************************
!
!!<subroutine>
!
!  subroutine optcana_analyticalError (rglobalData,rconstraints,rsolution,rreference,&
!      DerrorU,DerrorP,DerrorLambda,DerrorXi,boutput)
!
!!<description>
!  ! Computes the L2-error $||y-y0||_2$ of the given solution y to a reference
!  ! solution y0 given as analytical function.
!!</description>
!  
!!<input>
!  ! Solution vector to compute the norm/error from.
!  type(t_spacetimeVector), intent(IN) :: rsolution
!  
!  ! Analytic solution defining the reference function z.
!  type(t_anSolution), intent(inout) :: rreference
!
!  ! Constraints in the optimal control problem.
!  type(t_optcconstraintsSpaceTime), intent(in) :: rconstraints
!  
!  ! Global settings for callback routines.
!  type(t_globalData), intent(inout), target :: rglobalData
!  
!  ! Flag that determines if the error in each timestep is written to the terminal.
!  logical, intent(in) :: boutput
!!</input>
!
!!<output>
!  ! Returns information about the error. ||y-y0||_{L^2}.
!  ! DerrorU(1) = on the interval [0,T]
!  ! DerrorU(2) = on the interval [0,T)
!  ! DerrorU(3) = on the interval (0,T]
!  ! DerrorU(4) = on the interval (0,T)
!  real(DP), dimension(:), intent(out) :: DerrorU
!
!  ! Returns information about the error in the pressure. ||p-p0||_{L^2}.
!  ! DerrorP(1) = on the interval [0,T]
!  ! DerrorP(2) = on the interval [0,T)
!  ! DerrorP(3) = on the interval (0,T]
!  ! DerrorP(4) = on the interval (0,T)
!  real(DP), dimension(:), intent(out) :: DerrorP
!
!  ! Returns information about the error. ||lambda-lambda0||_{L^2}.
!  ! DerrorLambda(1) = on the interval [0,T]
!  ! DerrorLambda(2) = on the interval [0,T)
!  ! DerrorLambda(3) = on the interval (0,T]
!  ! DerrorLambda(4) = on the interval (0,T)
!  real(DP), dimension(:), intent(out) :: DerrorLambda
!
!  ! Returns information about the error in the pressure. ||xi-xi0||_{L^2}.
!  ! DerrorXi(1) = on the interval [0,T]
!  ! DerrorXi(2) = on the interval [0,T)
!  ! DerrorXi(3) = on the interval (0,T]
!  ! DerrorXi(4) = on the interval (0,T)
!  real(DP), dimension(:), intent(out) :: DerrorXi
!!</output>
!  
!!</subroutine>
!    
!    ! local variables
!    integer :: isubstep,i
!    real(DP) :: dtstep,dtimePrimal,dtimeDual
!    real(DP) :: derrU, derrP, derrLambda, derrXi
!    real(DP),dimension(6) :: Derr
!    type(t_collection) :: rcollection
!    type(t_collection), target :: ruserCollection
!    type(t_vectorBlock) :: rtempVector
!    real(dp), dimension(:), pointer :: p_Ddata
!    
!    ! Create a temp vector
!    call lsysbl_createVectorBlock(rsolution%p_rspaceDiscr,rtempVector,.true.)
!    call lsysbl_getbase_double (rtempVector,p_Ddata)
!
!    ! Some basic initialisation.
!    ! Initialise the collection for the assembly process with callback routines.
!    ! This stores the simulation time in the collection and sets the
!    ! current subvector z for the callback routines.
!    call collct_init(rcollection)
!    
!    DerrorU(:) = 0.0_DP
!    DerrorP(:) = 0.0_DP
!    DerrorLambda(:) = 0.0_DP
!    DerrorXi(:) = 0.0_DP
!
!    do isubstep = 1,rsolution%NEQtime
!    
!      ! Current point in time.
!      ! For the first iterate, take a look at time interval 1 to
!      ! get the length and start point of the interval.
!      ! For all further iterates, look at the time interval leading to
!      ! that iterate.
!      if (isubstep .gt. 1) then
!        call tdiscr_getTimestep(rsolution%p_rtimeDiscr,isubstep-1,dtimePrimal,dtstep)
!        dtimeDual = dtimePrimal - (1.0_DP-rsolution%p_rtimeDiscr%dtheta)*dtstep
!      else
!        call tdiscr_getTimestep(rsolution%p_rtimeDiscr,isubstep,dtstep=dtstep,dtimestart=dtimePrimal)
!        dtimeDual = dtimePrimal
!      end if
!
!      ! Get the solution.
!      ! Evaluate the space time function in rvector in the point
!      ! in time dtime. Independent of the discretisation in time,
!      ! this will give us a vector in space.
!      ! Note 1: The dual solution is shifted by (1-dtheta)*dtstep!
!      ! We therefore only have to evaluate once!
!      ! Note 2: For Crank-Nicolson, the time discretisation discretises
!      ! the primal pressure at the point of the dual velocity and
!      ! the dual pressure at the time of the primal velocity.
!      ! We compensate for this time shift during the error calculation.
!      
!      !CALL sptivec_getTimestepData (rsolution, 1+isubstep, rtempVector)
!      call tmevl_evaluate(rsolution,dtimePrimal,rtempVector)
!
!      ! Use our standard implementation to evaluate the error.
!      if (rreference%ctype .ne. ANSOL_TP_ANALYTICAL) then
!        call ansol_prepareEval (rreference,rcollection,"SOL",dtimePrimal)
!      else
!        ! Prepare the user-defined collection for assembly.
!        call collct_init(ruserCollection)
!      end if
!  
!      ! Save the function type to the collection, so the callback knows how
!      ! to evaluate.
!      rcollection%IquickAccess(2) = rreference%ctype
!      
!      ! The user-defined collection is the follower of rcollection.
!      rcollection%p_rnextCollection => ruserCollection
!
!      ! Perform error analysis to calculate and add 1/2||y-y0||^2_{L^2},...
!      ! Primal velocity, dual pressure
!      do i=1,2
!        rcollection%IquickAccess(1) = i
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          ! Analytically given data
!          call user_initCollectForVecAssembly (rglobalData,&
!              rreference%iid,i,dtimePrimal,ruserCollection)
!        end if
!
!        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
!            optcana_evalFunction,rcollection)
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          call user_doneCollectForVecAssembly (rglobalData,ruserCollection)
!        end if
!      end do
!      
!      ! Primal pressure only in the 1st timestep.
!      if (isubstep .eq. 0) then
!        i=3
!        rcollection%IquickAccess(1) = i
!        
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          ! Analytically given data
!          call user_initCollectForVecAssembly (rglobalData,&
!              rreference%iid,i,dtimePrimal,ruserCollection)
!        end if
!        
!        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
!            optcana_evalFunction,rcollection)
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          call user_doneCollectForVecAssembly (rglobalData,ruserCollection)
!        end if
!      end if
!      
!      if (isubstep .ne. rsolution%NEQtime) then
!        i=6
!        rcollection%IquickAccess(1) = i
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          ! Analytically given data
!          call user_initCollectForVecAssembly (rglobalData,&
!              rreference%iid,i,dtimePrimal,ruserCollection)
!        end if
!
!        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
!            optcana_evalFunction,rcollection)
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          call user_doneCollectForVecAssembly (rglobalData,ruserCollection)
!        end if
!      end if
!          
!      ! The same for the dual equation.
!      ! Dual velocity, primal pressure.
!      ! In rtempVector(4..6) is the dual solution at time dtimeDual,
!      ! so we don"t have to evaluate the function again!
!
!      ! If we have a function in rreference, switch the time for ir.
!      if (rreference%ctype .ne. ANSOL_TP_ANALYTICAL) then
!
!        call ansol_doneEval (rcollection,"SOL")
!        call ansol_prepareEval (rreference,rcollection,"SOL",dtimeDual)
!
!      end if
!
!      do i=4,5
!        rcollection%IquickAccess(1) = i
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          ! Analytically given data
!          call user_initCollectForVecAssembly (rglobalData,&
!              rreference%iid,i,dtimeDual,ruserCollection)
!        end if
!
!        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
!            optcana_evalFunction,rcollection)
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          call user_doneCollectForVecAssembly (rglobalData,ruserCollection)
!        end if
!            
!      end do
!      
!      if (isubstep .ne. 0) then
!        i=3
!        rcollection%IquickAccess(1) = i
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          ! Analytically given data
!          call user_initCollectForVecAssembly (rglobalData,&
!              rreference%iid,i,dtimeDual,ruserCollection)
!        end if
!
!        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
!            optcana_evalFunction,rcollection)
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          call user_doneCollectForVecAssembly (rglobalData,ruserCollection)
!        end if
!      end if
!          
!      ! Dual pressure only in the last timestep.
!      if (isubstep .eq. rsolution%NEQtime) then
!        i=6
!        rcollection%IquickAccess(1) = i
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          ! Analytically given data
!          call user_initCollectForVecAssembly (rglobalData,&
!              rreference%iid,i,dtimeDual,ruserCollection)
!        end if
!
!        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
!            optcana_evalFunction,rcollection)
!
!        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!          call user_doneCollectForVecAssembly (rglobalData,ruserCollection)
!        end if
!      end if
!
!      ! Clean up the collection -- either ours or the user-defined one.
!      if (rreference%ctype .ne. ANSOL_TP_ANALYTICAL) then
!        call ansol_doneEval (rcollection,"SOL")
!      else
!        call collct_done (ruserCollection)
!      end if
!
!      ! Get the errors in that timestep.
!      derrU = sqrt(Derr(1)**2 + Derr(2)**2)
!      derrP = Derr(3)
!      derrLambda = sqrt(Derr(4)**2 + Derr(5)**2)
!      derrXi = Derr(6)
!
!      ! We use the summed trapezoidal rule.
!      ! Watch out with the start/end of the time interval when
!      ! plugging the values into the error arrays.
!      if (isubstep .eq. 1) then
!        
!        DerrorU(1)      = DerrorU(1)      + 0.5_DP*derrU**2
!        DerrorP(1)      = DerrorP(1)      + 0.5_DP*derrP**2
!        DerrorLambda(1) = DerrorLambda(1) + 0.5_DP*derrLambda**2
!        DerrorXi(1)     = DerrorXi(1)     + 0.5_DP*derrXi**2
!
!        DerrorU(2)      = DerrorU(2)      + 0.5_DP*derrU**2
!        DerrorP(2)      = DerrorP(2)      + 0.5_DP*derrP**2
!        DerrorLambda(2) = DerrorLambda(2) + 0.5_DP*derrLambda**2
!        DerrorXi(2)     = DerrorXi(2)     + 0.5_DP*derrXi**2
!      
!      end if
!
!      if (isubstep .eq. 2) then
!        
!        DerrorU(1)      = DerrorU(1)      + derrU**2
!        DerrorP(1)      = DerrorP(1)      + derrP**2
!        DerrorLambda(1) = DerrorLambda(1) + derrLambda**2
!        DerrorXi(1)     = DerrorXi(1)     + derrXi**2
!
!        DerrorU(2)      = DerrorU(2)      + derrU**2
!        DerrorP(2)      = DerrorP(2)      + derrP**2
!        DerrorLambda(2) = DerrorLambda(2) + derrLambda**2
!        DerrorXi(2)     = DerrorXi(2)     + derrXi**2
!      
!        DerrorU(3)      = DerrorU(3)      + 0.5_DP*derrU**2
!        DerrorP(3)      = DerrorP(3)      + 0.5_DP*derrP**2
!        DerrorLambda(3) = DerrorLambda(3) + 0.5_DP*derrLambda**2
!        DerrorXi(3)     = DerrorXi(3)     + 0.5_DP*derrXi**2
!
!        DerrorU(4)      = DerrorU(4)      + 0.5_DP*derrU**2
!        DerrorP(4)      = DerrorP(4)      + 0.5_DP*derrP**2
!        DerrorLambda(4) = DerrorLambda(4) + 0.5_DP*derrLambda**2
!        DerrorXi(4)     = DerrorXi(4)     + 0.5_DP*derrXi**2
!      
!      end if
!
!      if ((isubstep .ge. 3) .and. (isubstep .le. rsolution%NEQtime-2)) then
!        
!        DerrorU(1)      = DerrorU(1)      + derrU**2
!        DerrorP(1)      = DerrorP(1)      + derrP**2
!        DerrorLambda(1) = DerrorLambda(1) + derrLambda**2
!        DerrorXi(1)     = DerrorXi(1)     + derrXi**2
!
!        DerrorU(2)      = DerrorU(2)      + derrU**2
!        DerrorP(2)      = DerrorP(2)      + derrP**2
!        DerrorLambda(2) = DerrorLambda(2) + derrLambda**2
!        DerrorXi(2)     = DerrorXi(2)     + derrXi**2
!      
!        DerrorU(3)      = DerrorU(3)      + derrU**2
!        DerrorP(3)      = DerrorP(3)      + derrP**2
!        DerrorLambda(3) = DerrorLambda(3) + derrLambda**2
!        DerrorXi(3)     = DerrorXi(3)     + derrXi**2
!
!        DerrorU(4)      = Derroru(4)      + derrU**2
!        DerrorP(4)      = DerrorP(4)      + derrP**2
!        DerrorLambda(4) = DerrorLambda(4) + derrLambda**2
!        DerrorXi(4)     = DerrorXi(4)     + derrXi**2
!      
!      end if
!      
!      if (isubstep .eq. rsolution%NEQtime-1) then
!        
!        DerrorU(1)      = DerrorU(1)      + derrU**2
!        DerrorP(1)      = DerrorP(1)      + derrP**2
!        DerrorLambda(1) = DerrorLambda(1) + derrLambda**2
!        DerrorXi(1)     = DerrorXi(1)     + derrXi**2
!
!        DerrorU(2)      = DerrorU(2)      + 0.5_DP*derrU**2
!        DerrorP(2)      = DerrorP(2)      + 0.5_DP*derrP**2
!        DerrorLambda(2) = DerrorLambda(2) + 0.5_DP*derrLambda**2
!        DerrorXi(2)     = DerrorXi(2)     + 0.5_DP*derrXi**2
!      
!        DerrorU(3)      = DerrorU(3)      + derrU**2
!        DerrorP(3)      = DerrorP(3)      + derrP**2
!        DerrorLambda(3) = DerrorLambda(3) + derrLambda**2
!        DerrorXi(3)     = DerrorXi(3)     + derrXi**2
!
!        DerrorU(4)      = DerrorU(4)      + 0.5_DP*derrU**2
!        DerrorP(4)      = DerrorP(4)      + 0.5_DP*derrP**2
!        DerrorLambda(4) = DerrorLambda(4) + 0.5_DP*derrLambda**2
!        DerrorXi(4)     = DerrorXi(4)     + 0.5_DP*derrXi**2
!      
!      end if
!
!      if (isubstep .eq. rsolution%NEQtime) then
!        
!        DerrorU(1)      = DerrorU(1)      + 0.5_DP*derrU**2
!        DerrorP(1)      = DerrorP(1)      + 0.5_DP*derrP**2
!        DerrorLambda(1) = DerrorLambda(1) + 0.5_DP*derrLambda**2
!        DerrorXi(1)     = DerrorXi(1)     + 0.5_DP*derrXi**2
!
!        DerrorU(3)      = DerrorU(3)      + 0.5_DP*derrU**2
!        DerrorP(3)      = DerrorP(3)      + 0.5_DP*derrP**2
!        DerrorLambda(3) = DerrorLambda(3) + 0.5_DP*derrLambda**2
!        DerrorXi(3)     = DerrorXi(3)     + 0.5_DP*derrXi**2
!
!      end if
!      
!      if (boutput) then
!        call output_line("error("//trim(sys_siL(isubstep,10))//") = "// &
!            trim(sys_sdEL(Derr(1),10))//" / "//&
!            trim(sys_sdEL(Derr(2),10))//" / "// &
!            trim(sys_sdEL(Derr(3),10))//" / "// &
!            trim(sys_sdEL(Derr(4),10))//" / "//&
!            trim(sys_sdEL(Derr(5),10))//" / "// &
!            trim(sys_sdEL(Derr(6),10)))
!      end if
!          
!    end do
!
!    ! Get the error return values.
!    DerrorU = sqrt(DerrorU*dtstep)
!    DerrorP = sqrt(DerrorP*dtstep)
!    DerrorLambda = sqrt(DerrorLambda*dtstep)
!    DerrorXi = sqrt(DerrorXi*dtstep)
!
!    ! Release temnp vector and temp data
!    call lsysbl_releaseVector (rtempVector)
!    call collct_done(rcollection)
!
!  end subroutine

!hier weiter

end module
