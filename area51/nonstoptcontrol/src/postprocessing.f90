!##############################################################################
!# ****************************************************************************
!# <name> spatialoperators </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains routines for assembling and applying operators in space.
!# </purpose>
!##############################################################################

module postprocessing

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  use derivatives
  use triangulation
  use pprocerror
  
  use ucd
 
  use physics
  use spacetimevectors
  use spacetimematvec
  
  implicit none

contains

!<subroutine>

  subroutine ferrFunction (cderivative, rdiscretisation, &
                                  nelements, npointsPerElement, Dpoints, &
                                  IdofsTest, rdomainIntSubset, &
                                  Dvalues, rcollection)
  
  use fsystem
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
!<description>
  ! Callback function to calculate the error of a FE function to the
  ! corresponding analytical function.
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
    integer :: icomponent,cequation,creferenceProblem,ierrType
    real(DP) :: dtime,dalpha

    ! Get the settings of the equation
    icomponent = rcollection%IquickAccess (1)
    cequation = rcollection%IquickAccess (2)
    creferenceProblem = rcollection%IquickAccess (3)
    dtime = rcollection%DquickAccess (1)
    dalpha = rcollection%DquickAccess (2)

    if (cderivative .eq. DER_FUNC) then
      ! Select the equation and calculate the error.
      select case (cequation)
      case (0,2)
        ! -------------------------------------------------------------
        ! Heat equation
        ! -------------------------------------------------------------
        select case (icomponent)
        case (1)
          !Dvalues(:,:) = dtime*(Dpoints(1,:,:)*Dpoints(2,:,:))
          !Dvalues(:,:) = 0.0_DP
          
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_heatY1 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
          case (2)
            ! 2.) -> BC in spacetimebc.f90 beachten!
            Dvalues(:,:) = fct_heatY2 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (3)
            ! 3.) -> BC in spacetimebc.f90 beachten!
            Dvalues(:,:) = fct_heatY3 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (4)
            ! 4.) -> BC in spacetimebc.f90 beachten!
            Dvalues(:,:) = fct_heatY4 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (5)
            ! 5.) -> BC in spacetimebc.f90 beachten!
            Dvalues(:,:) = fct_heatY5 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (6)
            ! 6.) -> BC in spacetimebc.f90 beachten!
            Dvalues(:,:) = fct_heatY6 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
          
        case (2)
          !Dvalues(:,:) = (1.0_DP-dtime)*(Dpoints(1,:,:)*Dpoints(2,:,:))
          !Dvalues(:,:) = 0.0_DP
          !Dvalues(:,:) = - (2.0_DP*dtime**2 * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) 
          ! +  2.0_DP*dtime**2.0_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )
          
          ! For 1D: See BC in spacetimebc.f90!!!
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_heatLambda1 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
            
          case (2)
            ! 2.) 
            Dvalues(:,:) = fct_heatLambda2 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (3)
            ! 3.) 
            Dvalues(:,:) = fct_heatLambda3 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
            
          case (4)
            ! 4.) 
            Dvalues(:,:) = fct_heatLambda4 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
            
          case (5)
            ! 5.) 
            Dvalues(:,:) = fct_heatLambda5 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
            
          case (6)
            ! 6.) 
            Dvalues(:,:) = fct_heatLambda6 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
            
          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
        case default
          ! Should not happen
          call sys_halt()
        end select

      case (1)
        ! -------------------------------------------------------------
        ! Stokes equation
        ! -------------------------------------------------------------
        select case (icomponent)
        
        ! Primal BC
        case (1)
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesY1_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesY2_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesY3_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesY4_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesY5_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesY6_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesY7_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
          
        case (2)
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesY1_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesY2_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesY3_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesY4_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesY5_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesY6_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesY7_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
        
        case (3)
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesP1 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesP2 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesP3 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesP4 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesP5 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesP6 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesP7 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
        
        ! Dual BC
        case (4)
        
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesLambda1_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesLambda2_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesLambda3_x (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesLambda4_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesLambda5_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesLambda6_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesLambda7_x(Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
          
        case (5)
        
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesLambda1_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesLambda2_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)
          
          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesLambda3_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesLambda4_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesLambda5_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesLambda6_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesLambda7_y (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
        
        case (6)
        
          select case (creferenceProblem)
          case (0)
          case (1)
            ! 1.)
            Dvalues(:,:) = fct_stokesXi1 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (2)
            ! 2.)
            Dvalues(:,:) = fct_stokesXi2 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (3)
            ! 3.)
            Dvalues(:,:) = fct_stokesXi3 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (4)
            ! 4.)
            Dvalues(:,:) = fct_stokesXi4 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (5)
            ! 5.)
            Dvalues(:,:) = fct_stokesXi5 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (6)
            ! 6.)
            Dvalues(:,:) = fct_stokesXi6 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case (7)
            ! 7.)
            Dvalues(:,:) = fct_stokesXi7 (Dpoints(1,:,:),Dpoints(2,:,:),dtime,dalpha)

          case default
            call output_line ("Problem not supported.")
            call sys_halt()
          end select
        
        case default
          ! Should not happen
          call sys_halt()
        end select
        
      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()
        
      end select

    else
    
      call output_line ("only function values supported.")
      call sys_halt()
      
    end if

  end subroutine 
  
  ! ***************************************************************************

  subroutine stpp_postproc (rphysics,rvector,bwriteUCD,bcalcError)

    ! Postprocessing of a space-time vector.
    
    ! Structure defining the problem physics
    type(t_physics), intent(in) :: rphysics
    
    ! Vector to be postprocessed.
    type(t_spacetimeVector), intent(in) :: rvector
    
    ! Whether to write UCD output files or not.
    logical, intent(in) :: bwriteUCD

    ! Whether to calculate the error to an analytical reference function or not.
    logical, intent(in) :: bcalcError
  
    ! local variables
    integer :: istep
    type(t_vectorBlock) :: rspaceVec
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    type(t_triangulation), pointer :: p_rtria
    real(DP), dimension(6) :: Derrors
    real(DP), dimension(6) :: DerrorTotal
    type(t_collection) :: rcollection
    real(DP) :: dtimePrimal,dtimeDual,dtstep
    integer :: i

    call lsysbl_createVectorblock (rvector%p_rspaceDiscr,rspaceVec)

    p_rtria => rspaceVec%p_rblockDiscr%p_rtriangulation

    DerrorTotal(:) = 0.0_DP
    Derrors(:) = 0.0_DP

    do istep=1,rvector%NEQtime
      
      ! CN -> Dual is at different time!
      call tdiscr_getTimestep(rvector%p_rtimeDiscr,istep-1,dtimePrimal,dtstep)
      if (istep .eq. 1) then
        dtimeDual = dtimePrimal
      else
        dtimeDual = dtimePrimal - (1.0_DP-rvector%p_rtimeDiscr%dtheta)*dtstep
      end if
      
      select case (rphysics%cequation)
      case (0,2)
        ! Heat equation
      
        call sptivec_getTimestepData (rvector, istep, rspaceVec)

        if (bwriteUCD) then
          ! Write to vtk file.
          
          call ucd_startVTK(rexport,UCD_FLAG_STANDARD,&
              rvector%p_rspaceDiscr%p_rtriangulation,&
              "./gmv/u.vtk."//trim(sys_si0L(istep-1,5)))
              
          call lsyssc_getbase_double (rspaceVec%RvectorBlock(1),p_Ddata1)
          call ucd_addVariableVertexBased(rexport,"y",UCD_VAR_STANDARD,p_Ddata1(1:p_rtria%NVT))

          call lsyssc_getbase_double (rspaceVec%RvectorBlock(2),p_Ddata1)
          call ucd_addVariableVertexBased(rexport,"lambda",UCD_VAR_STANDARD,p_Ddata1(1:p_rtria%NVT))
          
          call ucd_write(rexport)
          call ucd_release (rexport)
        end if

        if (bcalcError) then
          rcollection%IquickAccess (2) = rphysics%cequation
          rcollection%IquickAccess (3) = rphysics%creferenceProblem
          rcollection%DquickAccess (2) = rphysics%doptControlAlpha

          rcollection%DquickAccess (1) = dtimePrimal
          rcollection%IquickAccess (1) = 1
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)),&
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%DquickAccess (1) = dtimeDual
          rcollection%IquickAccess (1) = 2
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)),&
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          call output_line ("Error("//trim(sys_siL(istep,10))//") = "// &
              trim(sys_sdEL(Derrors(1),5))//" "// &
              trim(sys_sdEL(Derrors(2),5)) )
              
          ! Summed trapezoidal rule, without first and last error.
          if ((istep .ne. 1) .and. (istep .ne. rvector%NEQtime)) then
            if ((istep .eq. 2) .or. (istep .eq. rvector%NEQtime-1)) then
              DerrorTotal(:) = DerrorTotal(:) + 0.5_DP * Derrors(:)**2
            else
              DerrorTotal(:) = DerrorTotal(:) + Derrors(:)**2
            end if
          end if
          
          if (istep .eq. rvector%NEQtime) then
            call output_lbrk()
            call output_line ("Error(Total) = "// &
                trim(sys_sdEL(sqrt(DerrorTotal(1)*dtstep),5))//" "// &
                trim(sys_sdEL(sqrt(DerrorTotal(2)*dtstep),5)) )
          end if
          
        end if
        
      case (1)
        ! Stokes equation
        
        call sptivec_getTimestepData (rvector, istep, rspaceVec)

        if (bwriteUCD) then
          ! Write to vtk file.
          
          call ucd_startVTK(rexport,UCD_FLAG_STANDARD,&
              rvector%p_rspaceDiscr%p_rtriangulation,&
              "./gmv/u.vtk."//trim(sys_si0L(istep-1,5)))
              
          call lsyssc_getbase_double (rspaceVec%RvectorBlock(1),p_Ddata1)
          call lsyssc_getbase_double (rspaceVec%RvectorBlock(2),p_Ddata2)

          call ucd_addVarVertBasedVec(rexport,"y",p_Ddata1(1:p_rtria%NVT),p_Ddata2(1:p_rtria%NVT))

          call lsyssc_getbase_double (rspaceVec%RvectorBlock(3),p_Ddata1)
          call ucd_addVariableElementBased(rexport,"p",UCD_VAR_STANDARD,p_Ddata1(1:p_rtria%NEL))

          call lsyssc_getbase_double (rspaceVec%RvectorBlock(4),p_Ddata1)
          call lsyssc_getbase_double (rspaceVec%RvectorBlock(5),p_Ddata2)
          call ucd_addVarVertBasedVec(rexport,"lambda",p_Ddata1(1:p_rtria%NVT),p_Ddata2(1:p_rtria%NVT))
          
          call lsyssc_getbase_double (rspaceVec%RvectorBlock(6),p_Ddata1)
          call ucd_addVariableElementBased(rexport,"xi",UCD_VAR_STANDARD,p_Ddata1(1:p_rtria%NEL))
          
          call ucd_write(rexport)
          call ucd_release (rexport)
        end if

        if (bcalcError) then
          rcollection%IquickAccess (2) = rphysics%cequation
          rcollection%IquickAccess (3) = rphysics%creferenceProblem
          rcollection%DquickAccess (2) = rphysics%doptControlAlpha

          rcollection%DquickAccess (1) = dtimePrimal
          rcollection%IquickAccess (1) = 1
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%IquickAccess (1) = 2
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          rcollection%IquickAccess (1) = 6
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%DquickAccess (1) = dtimeDual
          rcollection%IquickAccess (1) = 3
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%IquickAccess (1) = 4
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          rcollection%IquickAccess (1) = 5
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          call output_line ("Error("//trim(sys_siL(istep,10))//") = "// &
              trim(sys_sdEL(Derrors(1),5))//" "// &
              trim(sys_sdEL(Derrors(2),5))//" "// &
              trim(sys_sdEL(Derrors(3),5))//" "// &
              trim(sys_sdEL(Derrors(4),5))//" "// &
              trim(sys_sdEL(Derrors(5),5))//" "// &
              trim(sys_sdEL(Derrors(6),5)) )
              
          ! Summed trapezoidal rule, without first and last error.
          if ((istep .ne. 1) .and. (istep .ne. rvector%NEQtime)) then
            if ((istep .eq. 2) .or. (istep .eq. rvector%NEQtime-1)) then
              DerrorTotal(:) = DerrorTotal(:) + 0.5_DP * Derrors(:)**2
            else
              DerrorTotal(:) = DerrorTotal(:) + Derrors(:)**2
            end if
          end if
          
          if (istep .eq. rvector%NEQtime) then
            call output_lbrk()
            call output_line ("Error(Total) = "// &
                trim(sys_sdEL(sqrt(DerrorTotal(1)*dtstep),5))//" "// &
                trim(sys_sdEL(sqrt(DerrorTotal(2)*dtstep),5))//" "// &
                trim(sys_sdEL(sqrt(DerrorTotal(3)*dtstep),5))//" "// &
                trim(sys_sdEL(sqrt(DerrorTotal(4)*dtstep),5))//" "// &
                trim(sys_sdEL(sqrt(DerrorTotal(5)*dtstep),5))//" "// &
                trim(sys_sdEL(sqrt(DerrorTotal(6)*dtstep),5)) )
          end if
              
        end if

      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()
      
      end select
      
    end do
    
    call lsysbl_releaseVector (rspaceVec)
  
  end subroutine
  
  ! ***************************************************************************

  subroutine stpp_printDefectSubnorms (rmatrix,rx,rb,rd)

    ! Prints the norms of the defects in each timestep to the terminal.
    
    ! Underlying matrix
    type(t_spaceTimeMatrix), intent(in) :: rmatrix
    
    ! Solution vector
    type(t_spacetimeVector), intent(in) :: rx
    
    ! RHS vector
    type(t_spacetimeVector), intent(in) :: rb
  
    ! Temporary vector
    type(t_spacetimeVector), intent(inout) :: rd
  
    ! local variables
    integer :: istep, icomp
    real(DP) :: dnormsum
    real(DP), dimension(:), allocatable :: Dnorms
    integer, dimension(:), allocatable :: Cnorms
    type(t_vectorBlock) :: rspacetemp
    real(DP), dimension(:), pointer :: p_Dy
    
    ! create a temp vector
    call lsysbl_createVectorBlock (rmatrix%p_rspaceDiscr,rspacetemp)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rspacetemp,p_Dy)
    
    ! Create the defect.
    call sptivec_copyVector (rb,rd)
    call stmv_matvec (rmatrix, rx, rd, -1.0_DP, 1.0_DP)
    call spop_applyBC (rmatrix%p_rboundaryCond, SPOP_DEFECT, rd)
    
    ! Print the defect norm in each timestep
    allocate (Dnorms(rspacetemp%nblocks))
    allocate (Cnorms(rspacetemp%nblocks))
    Cnorms(:) = LINALG_NORMEUCLID
    
    do istep=1,rd%NEQtime

      ! Calculate the euclidian norm
      call sptivec_getTimestepData (rd, istep, rspacetemp)
      call lsysbl_vectorNormBlock (rspacetemp,Cnorms,Dnorms)
      
      ! Calculate the L2-norm(s) of every block.
      dnormsum = sum(Dnorms)
      call output_line ("||D_"//trim(sys_siL(istep,10))//"|| = "//&
        trim(sys_sdEP(dnormsum/sqrt(real(rspaceTemp%NEQ,DP)),16,8)))
      do icomp = 1,rspacetemp%nblocks
        call output_line ("  ||D_"//trim(sys_siL(istep,10))//&
           "^"//trim(sys_siL(icomp,10))//"|| = "//&
           trim(sys_sdEP(Dnorms(icomp)/ &
                sqrt(real(rspaceTemp%RvectorBlock(icomp)%NEQ,DP)),16,8)))
      end do
    end do
    
    deallocate(Dnorms)
    
    !vergleich mit dem optc-Code: Ab der 2. Defektkomponente stimmt der Defekt nicht mehr überein
    
    call lsysbl_releaseVector (rspacetemp)
    
  end subroutine

  ! ***************************************************************************

  subroutine stpp_printDefectSubnormsDirect (rd)

    ! Prints the norms of the defects in each timestep to the terminal.
    
    ! Defect vector
    type(t_spacetimeVector), intent(inout) :: rd
  
    ! local variables
    integer :: istep, icomp
    real(DP) :: dnormsum
    real(DP), dimension(:), allocatable :: Dnorms
    integer, dimension(:), allocatable :: Cnorms
    type(t_vectorBlock) :: rspacetemp
    
    ! create a temp vector
    call lsysbl_createVectorBlock (rd%p_rspaceDiscr,rspacetemp)
    
    ! Print the defect norm in each timestep
    allocate (Dnorms(rspacetemp%nblocks))
    allocate (Cnorms(rspacetemp%nblocks))
    Cnorms(:) = LINALG_NORMEUCLID
    
    do istep=1,rd%NEQtime

      ! Calculate the euclidian norm
      call sptivec_getTimestepData (rd, istep, rspacetemp)
      call lsysbl_vectorNormBlock (rspacetemp,Cnorms,Dnorms)
      
      ! Calculate the L2-norm(s) of every block.
      dnormsum = sum(Dnorms)
      call output_line ("||D_"//trim(sys_siL(istep,10))//"|| = "//&
        trim(sys_sdEP(dnormsum/real(rspaceTemp%NEQ,DP),16,8)))
      do icomp = 1,rspacetemp%nblocks
        call output_line ("  ||D_"//trim(sys_siL(istep,10))//&
           "^"//trim(sys_siL(icomp,10))//"|| = "//&
           trim(sys_sdEP(Dnorms(icomp)/ &
                real(rspaceTemp%RvectorBlock(icomp)%NEQ,DP),16,8)))
      end do
    end do
    
    deallocate(Dnorms)
    
    call lsysbl_releaseVector (rspacetemp)
    
  end subroutine

end module
