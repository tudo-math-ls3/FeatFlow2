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
    integer :: ithetaschemetype

    call lsysbl_createVectorblock (rvector%p_rspaceDiscr,rspaceVec)

    p_rtria => rspaceVec%p_rblockDiscr%p_rtriangulation
    
    ithetaschemetype = rvector%p_rtimeDiscr%itag

    DerrorTotal(:) = 0.0_DP
    Derrors(:) = 0.0_DP

    do istep=1,rvector%NEQtime
      
      ! CN -> Dual is at different time!
      call tdiscr_getTimestep(rvector%p_rtimeDiscr,istep-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal

      ! Modified time scheme?
      if (ithetaschemetype .eq. 1) then
        if (istep .ne. 1) then
          dtimeDual = dtimePrimal - (1.0_DP-rvector%p_rtimeDiscr%dtheta)*dtstep
        end if
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
          rcollection%DquickAccess (3) = rphysics%dpar

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
          rcollection%DquickAccess (3) = rphysics%dpar

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
    
    !vergleich mit dem optc-Code: Ab der 2. Defektkomponente stimmt der Defekt nicht mehr �berein
    
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
