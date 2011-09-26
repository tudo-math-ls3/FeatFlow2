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
  use mprimitives
  use cubature
  use spdiscprojection
  
  use ucd
 
  use physics
  use spacetimevectors
  use spacetimematvec
  
  implicit none

contains

  ! ***************************************************************************

  subroutine stpp_postproc (rphysics,rvector,bwriteUCD,bcalcError,ccubTimeError,sucdfilename)

    ! Postprocessing of a space-time vector.
    
    ! Structure defining the problem physics
    type(t_physics), intent(in) :: rphysics
    
    ! Vector to be postprocessed.
    type(t_spacetimeVector), intent(in) :: rvector
    
    ! Whether to write UCD output files or not.
    logical, intent(in) :: bwriteUCD

    ! Whether to calculate the error to an analytical reference function or not.
    logical, intent(in) :: bcalcError
    
    ! Cubature formula to use for the time error.
    integer(I32), intent(in) :: ccubTimeError

    ! Name/Path of the UCD output file(s)
    character(len=*), intent(in) :: sucdfilename
  
    ! local variables
    integer :: istep
    type(t_vectorBlock) :: rspaceVec
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
    type(t_triangulation), pointer :: p_rtria
    real(DP), dimension(6) :: Derrors
    real(DP), dimension(6) :: DerrorTotal
    type(t_collection) :: rcollection
    real(DP) :: dtstep,dtimeY,dtimeLambda,dtimeP,dtimeXI
    integer :: ithetaschemetype

    call lsysbl_createVectorblock (rvector%p_rspaceDiscr,rspaceVec)

    p_rtria => rspaceVec%p_rblockDiscr%p_rtriangulation
    
    ithetaschemetype = rvector%p_rtimeDiscr%itag

    DerrorTotal(:) = 0.0_DP
    Derrors(:) = 0.0_DP

    do istep=1,rvector%NEQtime
      
      ! CN -> Dual is at different time!
      call tdiscr_getTimestep(rvector%p_rtimeDiscr,istep-1,dtimeY,dtstep)
      dtimeLambda = dtimeY

      ! Modified time scheme?
      select case (ithetaschemetype)
      case (0)
        ! Lambda at the same time as Y
        !
        ! Pressure p is half a  step before
        dtimeP = dtimeY
        if (istep .ne. 1) then
          dtimeP = dtimeY - (1.0_DP-rvector%p_rtimeDiscr%dtheta)*dtstep
        end if

        ! Pressure xi is half a step behind
        dtimeXI = dtimeY
        if (istep .ne. 1) then
          dtimeXI = dtimeY + (1.0_DP-rvector%p_rtimeDiscr%dtheta)*dtstep
        end if
        
      case (1)
        if (istep .ne. 1) then
          ! Lambda inbetween the Y
          dtimeLambda = dtimeY - (1.0_DP-rvector%p_rtimeDiscr%dtheta)*dtstep
        end if
        
        ! If the pressure is involved:
        dtimeP = dtimeLambda
        dtimeXI = dtimeY
      end select
      
      select case (rphysics%cequation)
      case (0,2)
        ! Heat equation
      
        call sptivec_getTimestepData (rvector, istep, rspaceVec)

        if (bwriteUCD) then
          ! Write to vtk file.
          
          call ucd_startVTK(rexport,UCD_FLAG_STANDARD,&
              rvector%p_rspaceDiscr%p_rtriangulation,&
              trim(sucdfilename)//"."//trim(sys_si0L(istep-1,5)))
              
          allocate (p_Ddata1(rvector%p_rspaceDiscr%p_rtriangulation%NVT))
              
          call spdp_projectToVertices (rspaceVec%RvectorBlock(1),p_Ddata1)
          call ucd_addVariableVertexBased(rexport,"y",UCD_VAR_STANDARD,p_Ddata1(1:p_rtria%NVT))

          call spdp_projectToVertices (rspaceVec%RvectorBlock(2),p_Ddata1)
          call ucd_addVariableVertexBased(rexport,"lambda",UCD_VAR_STANDARD,p_Ddata1(1:p_rtria%NVT))
          
          deallocate (p_Ddata1)
          
          call ucd_write(rexport)
          call ucd_release (rexport)
        end if

        if (bcalcError) then
          rcollection%IquickAccess (2) = rphysics%cequation
          rcollection%IquickAccess (3) = rphysics%creferenceProblem
          rcollection%DquickAccess (2) = rphysics%doptControlAlpha
          rcollection%DquickAccess (3) = rphysics%dpar
          rcollection%DquickAccess (4) = rphysics%dtimemin
          rcollection%DquickAccess (5) = rphysics%dtimemax

          rcollection%DquickAccess (1) = dtimeY
          rcollection%IquickAccess (1) = 1
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)),&
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%DquickAccess (1) = dtimeLambda
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
              trim(sucdfilename)//"."//trim(sys_si0L(istep-1,5)))
              
          allocate (p_Ddata1(rvector%p_rspaceDiscr%p_rtriangulation%NVT))
          allocate (p_Ddata2(rvector%p_rspaceDiscr%p_rtriangulation%NVT))
          allocate (p_Ddata3(rvector%p_rspaceDiscr%p_rtriangulation%NEL))
          
          call spdp_projectToVertices (rspaceVec%RvectorBlock(1),p_Ddata1)
          call spdp_projectToVertices (rspaceVec%RvectorBlock(2),p_Ddata2)

          call ucd_addVarVertBasedVec(rexport,"y",p_Ddata1(1:p_rtria%NVT),p_Ddata2(1:p_rtria%NVT))

          call spdp_projectToCells (rspaceVec%RvectorBlock(3),p_Ddata3)
          call ucd_addVariableElementBased(rexport,"p",UCD_VAR_STANDARD,p_Ddata3(1:p_rtria%NEL))

          call spdp_projectToVertices (rspaceVec%RvectorBlock(4),p_Ddata1)
          call spdp_projectToVertices (rspaceVec%RvectorBlock(5),p_Ddata2)
          call ucd_addVarVertBasedVec(rexport,"lambda",p_Ddata1(1:p_rtria%NVT),p_Ddata2(1:p_rtria%NVT))
          
          call spdp_projectToCells (rspaceVec%RvectorBlock(6),p_Ddata3)
          call ucd_addVariableElementBased(rexport,"xi",UCD_VAR_STANDARD,p_Ddata3(1:p_rtria%NEL))
          
          deallocate (p_Ddata3)
          deallocate (p_Ddata2)
          deallocate (p_Ddata1)
          
          call ucd_write(rexport)
          call ucd_release (rexport)
        end if

        if (bcalcError) then
          rcollection%IquickAccess (2) = rphysics%cequation
          rcollection%IquickAccess (3) = rphysics%creferenceProblem
          rcollection%DquickAccess (2) = rphysics%doptControlAlpha
          rcollection%DquickAccess (3) = rphysics%dpar
          rcollection%DquickAccess (4) = rphysics%dtimemin
          rcollection%DquickAccess (5) = rphysics%dtimemax

          rcollection%DquickAccess (1) = dtimeY
          rcollection%IquickAccess (1) = 1
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%IquickAccess (1) = 2
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          rcollection%DquickAccess (1) = dtimeP
          rcollection%IquickAccess (1) = 3
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%DquickAccess (1) = dtimeLambda
          rcollection%IquickAccess (1) = 4
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          rcollection%IquickAccess (1) = 5
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)), &
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          rcollection%DquickAccess (1) = dtimeXi
          rcollection%IquickAccess (1) = 6
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
    
    if (bcalcError) then
      call stpp_printError (rphysics,rvector,ccubTimeError)
    end if
  
  end subroutine
  
  ! ***************************************************************************

  subroutine stpp_printError (rphysics,rvector,ccubType)

    ! Prints the space-time error according to a specific 1D cubature formula.
    
    ! Structure defining the problem physics
    type(t_physics), intent(in) :: rphysics
    
    ! Vector to be postprocessed.
    type(t_spacetimeVector), intent(in) :: rvector
    
    ! Id of a 1D cubature formula which is used to calculate the error.
    integer(I32), intent(in) :: ccubType
  
    ! local variables
    integer :: istep
    type(t_vectorBlock) :: rspaceVec
    type(t_vectorBlock), pointer :: p_rx1,p_rx2
    integer :: ithetaschemetype
    real(DP) :: dtimePrimalStart,dtimePrimalEnd, dtimeDualStart,dtimeDualEnd,dtstep
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: ncubp,icubp
    real(DP), dimension(:,:), allocatable :: Dpoints
    real(DP), dimension(:), allocatable :: Domega
    real(DP), dimension(6) :: Derrors, DerrorsInterval
    real(DP), dimension(6) :: DerrorTotal
    real(DP) :: dpar
    type(t_collection) :: rcollection

    ! Create a temp vector pool for accessing space-time vectors.
    call sptivec_createAccessPool (rvector,raccessPool,2)
    call lsysbl_createVectorblock (rvector%p_rspaceDiscr,rspaceVec)

    ! Prepeate the cubature formula.
    ncubp = cub_igetNumPts(ccubType)
    allocate (Dpoints(NDIM1D,ncubp))
    allocate (Domega(ncubp))
    call cub_getCubature(ccubType, Dpoints, Domega)

    ithetaschemetype = rvector%p_rtimeDiscr%itag

    DerrorTotal(:) = 0.0_DP
    
    call output_lbrk()
    call output_line ("Error using cubature formula "//trim(sys_siL(ccubType,10))//":")

    do istep=1,rvector%NEQtime-1
      
      ! Get the interval points in time for the primal and dual equation;
      ! they may be different in case of CN!
      call tdiscr_getTimestep(rvector%p_rtimeDiscr,istep,&
          dtimePrimalEnd,dtstep,dtimePrimalStart)
      dtimeDualStart = dtimePrimalStart
      dtimeDualEnd = dtimePrimalEnd

      ! Modified time scheme?
      if (ithetaschemetype .eq. 1) then
        if (istep .ne. 1) then
          dtimeDualStart = dtimePrimalStart - (1.0_DP-rvector%p_rtimeDiscr%dtheta)*dtstep
        end if
        dtimeDualEnd = dtimePrimalEnd - (1.0_DP-rvector%p_rtimeDiscr%dtheta)*dtstep
      end if
      
      ! Read the timesteps at the ends of the interval.
      call sptivec_getVectorFromPool (raccessPool,istep,p_rx1)
      call sptivec_getVectorFromPool (raccessPool,istep+1,p_rx2)
      
      ! Perform a linear interpolation of the two vectors in time according to the
      ! cubature formula.
  
      DerrorsInterval(:) = 0.0_DP
  
      do icubp = 1,ncubp
        
        call mprim_linearRescale (Dpoints(1,icubp),-1.0_DP,1.0_DP,0.0_DP,1.0_DP,dpar)
        call lsysbl_copyVector (p_rx2,rspaceVec)
        call lsysbl_vectorLinearComb (p_rx1,rspaceVec,1.0_DP-dpar,dpar)
        
        ! Calculate the error in that point
        
        select case (rphysics%cequation)
        case (0,2)
          ! Heat equation
        
          rcollection%IquickAccess (2) = rphysics%cequation
          rcollection%IquickAccess (3) = rphysics%creferenceProblem
          rcollection%DquickAccess (2) = rphysics%doptControlAlpha
          rcollection%DquickAccess (3) = rphysics%dpar
          rcollection%DquickAccess (4) = rphysics%dtimemin
          rcollection%DquickAccess (5) = rphysics%dtimemax

          rcollection%DquickAccess (1) = (1.0_DP-dpar)*dtimePrimalStart + dpar*dtimePrimalEnd
          rcollection%IquickAccess (1) = 1
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)),&
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)

          rcollection%DquickAccess (1) = (1.0_DP-dpar)*dtimeDualStart + dpar*dtimeDualEnd
          rcollection%IquickAccess (1) = 2
          call pperr_scalar (PPERR_L2ERROR, Derrors(rcollection%IquickAccess (1)),&
              rspaceVec%RvectorBlock(rcollection%IquickAccess (1)),ferrFunction, &
              rcollection)
          
          call output_line ("  Error("//trim(sys_siL(istep,10))//"/"// &
              trim(sys_siL(icubp,10))//") = "// &
              trim(sys_sdEL(Derrors(1),5))//" "// &
              trim(sys_sdEL(Derrors(2),5)) )
              
        case (1)
          ! Stokes equation
          rcollection%IquickAccess (2) = rphysics%cequation
          rcollection%IquickAccess (3) = rphysics%creferenceProblem
          rcollection%DquickAccess (2) = rphysics%doptControlAlpha
          rcollection%DquickAccess (3) = rphysics%dpar
          rcollection%DquickAccess (4) = rphysics%dtimemin
          rcollection%DquickAccess (5) = rphysics%dtimemax
          
          rcollection%DquickAccess (1) = (1.0_DP-dpar)*dtimePrimalStart + dpar*dtimePrimalEnd
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

          rcollection%DquickAccess (1) = (1.0_DP-dpar)*dtimeDualStart + dpar*dtimeDualEnd
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
          
          call output_line ("  Error("//trim(sys_siL(istep,10))//"/"// &
              trim(sys_siL(icubp,10))//") = "// &
              trim(sys_sdEL(Derrors(1),5))//" "// &
              trim(sys_sdEL(Derrors(2),5))//" "// &
              trim(sys_sdEL(Derrors(3),5))//" "// &
              trim(sys_sdEL(Derrors(4),5))//" "// &
              trim(sys_sdEL(Derrors(5),5))//" "// &
              trim(sys_sdEL(Derrors(6),5)) )

        end select
      
        ! Sum up the error to the global error according to the cubature formula.
        DerrorsInterval(:) = DerrorsInterval(:) + &
            dtstep*0.5_DP*Domega(icubp) * Derrors(:)**2
            
      end do

      ! Sum up the interval error to the global error.
      if ((istep .ne. 1) .and. (istep .ne. rvector%NEQtime-1)) then
        DerrorTotal(:) = DerrorTotal(:) + DerrorsInterval(:)
      end if
      
      ! To create the norm of the error on the interval, take the square root.
      DerrorsInterval(:) = sqrt(DerrorsInterval(:))
      
      ! Print the interval error
      select case (rphysics%cequation)
      case (0,2)
        call output_line ("Error("//trim(sys_siL(istep,10))//") = "// &
            trim(sys_sdEL(DerrorsInterval(1),5))//" "// &
            trim(sys_sdEL(DerrorsInterval(2),5)) )
      case (1)
        call output_line ("Error("//trim(sys_siL(istep,10))//") = "// &
            trim(sys_sdEL(DerrorsInterval(1),5))//" "// &
            trim(sys_sdEL(DerrorsInterval(2),5))//" "// &
            trim(sys_sdEL(DerrorsInterval(3),5))//" "// &
            trim(sys_sdEL(DerrorsInterval(4),5))//" "// &
            trim(sys_sdEL(DerrorsInterval(5),5))//" "// &
            trim(sys_sdEL(DerrorsInterval(6),5)) )
      end select
      
    end do
    
    call lsysbl_releaseVector (rspaceVec)
    call sptivec_releaseAccessPool (raccessPool)
    
    deallocate (Domega)
    deallocate (Dpoints)
    
    ! To create the norm of the error, take the square root.
    DerrorTotal(:) = sqrt(DerrorTotal(:))
    
    ! Print the Total error
    call output_lbrk()
    select case (rphysics%cequation)
    case (0,2)
      call output_line ("Error(Total,cub) = "// &
          trim(sys_sdEL(DerrorTotal(1),5))//" "// &
          trim(sys_sdEL(DerrorTotal(2),5)) )
    case (1)
      call output_line ("Error(Total,cub) = "// &
          trim(sys_sdEL(DerrorTotal(1),5))//" "// &
          trim(sys_sdEL(DerrorTotal(2),5))//" "// &
          trim(sys_sdEL(DerrorTotal(3),5))//" "// &
          trim(sys_sdEL(DerrorTotal(4),5))//" "// &
          trim(sys_sdEL(DerrorTotal(5),5))//" "// &
          trim(sys_sdEL(DerrorTotal(6),5)) )
    end select
    call output_lbrk()
  
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
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    
    ! Get the discretiosation in space.
    call sth_getLevel(rmatrix%p_rspaceTimeHierarchy,rmatrix%ilevel,&
        p_rfeSpaceLevel=p_rfeSpaceLevel)
    
    ! create a temp vector
    call lsysbl_createVectorBlock (p_rfeSpaceLevel%p_rdiscretisation,rspacetemp)
    
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
