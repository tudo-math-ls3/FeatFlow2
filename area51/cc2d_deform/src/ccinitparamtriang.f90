!##############################################################################
!# ****************************************************************************
!# <name> ccinitparamtriang </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic initialisation routines for CC2D:
!# Initialisation of the parametrisation and the triangulations on all levels.
!#
!# 1.) cc_initParamTriang
!#     -> Read parametrisation, read triangulation, refine the mesh
!#
!# 2.) cc_doneParamTriang
!#     -> Remove the meshes of all levels from the heap
!#
!# </purpose>
!##############################################################################

module ccinitparamtriang

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use statistics
  use griddeform
  use collection
  use convection
  use ucd
    
  use ccbasic
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initParamTriang (rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the NLMIN/NLMAX parameters
  ! from the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ilvmin,ilvmax
  
    ! Variable for a filename:
    character(LEN=SYS_STRLEN) :: sString
    character(LEN=SYS_STRLEN) :: sPRMFile, sTRIFile
    type(t_timer) :: rtimer

    integer(i32), dimension(:), pointer :: p_InodalProperty
    integer(i32), dimension(:), pointer :: p_InodalProperty1

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)
    
    ! Get the .prm and the .tri file from the parameter list.
    ! note that parlst_getvalue_string returns us exactly what stands
    ! in the parameter file, so we have to apply READ to get rid of
    ! probable ''!
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sParametrisation',sString)
    read (sString,*) sPRMFile
                              
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sMesh',sString)
    read (sString,*) sTRIFile
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rproblem%rboundary, sPrmFile)
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation, &
        sTRIFile, rproblem%rboundary)

    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(rproblem%NLMIN-1,&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)

    ! Now, refine to level up to nlmax.
    do i=rproblem%NLMIN+1,rproblem%NLMAX
      call tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
    end do
    
    call storage_getbase_int (rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation%h_InodalProperty,&
    p_InodalProperty)
    
    call cc_deform(rproblem)
    
    ! Compress the level hierarchy.
    ! Share the vertex coordinates of all levels, so the coarse grid coordinates
    ! are 'contained' in the fine grid coordinates. The effect is:
    ! 1.) Save some memory
    ! 2.) Every change in the fine grid coordinates also affects the coarse
    !     grid coordinates and vice versa.
    do i=rproblem%NLMAX-1,rproblem%NLMIN,-1
      call tria_compress2LevelOrdHierarchy (rproblem%RlevelInfo(i+1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimeGridGeneration = &
      rproblem%rstatistics%dtimeGridGeneration + rtimer%delapsedReal

  end subroutine

! ***************************************************************************
  
!<subroutine>

  subroutine cc_deform(rproblem)
  
!<description>
  !
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

  ! Definitions of variables.
  !
  ! We need a couple of variables for this problem. Let's see...
  !
  ! An object for saving the domain:
  type(t_boundary) :: rboundary
  
  ! An object for saving the triangulation on the domain
  type(t_triangulation) :: rtriangulation

  ! An object specifying the discretisation.
  ! This contains also information about trial/test functions,...
  type(t_blockDiscretisation) :: rdiscretisation
  
  ! In this vector we will save the area of the deformed grid
  type(t_vectorScalar) :: rvectorAreaQ0

  ! NLMAX receives the level where we want to solve.
  integer :: NLMAX,NLMIN
  
  ! dummy
  integer :: idummy,iloop
  
  ! Output block for UCD output to GMV file
  type(t_ucdExport) :: rexport
  
  ! grid deform structure
  type(t_griddefInfo) :: rgriddefInfo

  ! Ok, let's start.
  !
  print *,"---------NumberOfElements---------"
  print *,rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation%NEL
  call output_lbrk ()


  ! griddeformation setup
  !rgriddefInfo,NLMIN,NLMAX,rboundary,iStyle,iadaptSteps
  call griddef_deformationInit(rgriddefInfo,&
                               rproblem%NLMIN,&
                               rproblem%NLMAX,rproblem%rboundary,GRIDDEF_CLASSICAL,5)
  
  do iloop=rproblem%NLMIN,rproblem%NLMAX
    call griddef_buildHGrid(rgriddefInfo,rproblem%RlevelInfo(iloop)%rtriangulation,iloop)
  end do
  
  ! Deform
  call griddef_performDeformation(rgriddefInfo,idummy,&
                                  .TRUE., .FALSE., .FALSE., &
                                  .FALSE., NLMAX, 0, 0,&
                                  tridef2d_monitorfct)
                                  
  ! Release Deformation structures
  call griddef_DeformationDone(rgriddefInfo)
    
  end subroutine

!******************************************************************************
  
  !<subroutine>
    subroutine tridef2d_monitorfct(DvertexCoords,Dentries)
  
  
    !<description>
      ! In this function we build the nodewise area distribution out
      ! of an elementwise distribution
    !</description>

    !<inputoutput>
     real(DP), dimension(:,:) :: DvertexCoords
     real(DP), dimension(:) :: Dentries
    !</inputoutput>

    !</subroutine>
    ! local variables
     real(dp),dimension(:,:),allocatable :: Dpoints
     integer :: ive,i1,ipoints
     integer :: iMethod
     real(DP) :: Dist,t,dt,dmin,cx,cy,rad1,rad2
     cx=0.5_dp
     cy=0.5_dp
     rad1=0.1_dp
     rad2=0.1_dp
     iMethod = 1
      
     ipoints = ubound(Dentries,1)
      
      
     select case(iMethod)
       ! case(0) 0.5 + x
       case(0)
         ! loop over all vertices and compute the monitor function
        do ive=1, ipoints
          Dentries(ive) = 0.5_dp + DvertexCoords(1,ive)
         end do
         
       ! Case(1) circle
       case(1)
         ! loop over all vertices and compute the monitor function
         do ive=1,ipoints
           Dist = sqrt((cx - DvertexCoords(1,ive))**2 + (cy - DvertexCoords(2,ive))**2)
           ! Good now define the monitor function
           Dist = abs(Dist - rad1)/rad1
           Dist=max(dist,0.04_dp)
           Dist=min(1.0_dp,dist)
           Dentries(ive)=Dist
         end do
       ! Case(2) do nothing
       case(2)
        do ive=1,ipoints
          Dentries(ive) = 1.0_dp
        end do
       ! Case(3) Elipse
       case(3)
       
        allocate(Dpoints(2,10000))
        dt = 6.28_dp/real(10000)
        t  = 0.0_dp
        do i1=1,10000
         Dpoints(1,i1) = cx + rad2 * cos(t)
         Dpoints(2,i1) = cy + rad1 * sin(t)
         t = t + dt
        end do
        ! loop over all points on the elipse
        do ive=1,ipoints
          dmin = 10000.0_dp
          do i1=1,10000
            Dist = sqrt((Dpoints(1,i1)-DvertexCoords(1,ive))**2 + (Dpoints(2,i1)-DvertexCoords(2,ive))**2)
            dmin =min(Dist,dmin)
          end do
          dmin=max(dmin,0.03_dp)
          dmin=min(1.0_dp,dmin)
          Dentries(ive) = dmin
        end do
       !Case 4
       case(4)
         do ive=1,ipoints
           Dist = sqrt((0.5_dp - DvertexCoords(1,ive))**2 + (0.2_dp - DvertexCoords(2,ive))**2)
           ! Good now define the monitor function
           Dist = abs(Dist - 0.05_dp)/0.05_dp
           Dist=max(dist,0.1_dp)
           Dist=min(1.0_dp,dist)
           Dentries(ive)=Dist
         end do
       case default
     end select
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    ! Release the triangulation on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      call tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)
    
  end subroutine

end module
