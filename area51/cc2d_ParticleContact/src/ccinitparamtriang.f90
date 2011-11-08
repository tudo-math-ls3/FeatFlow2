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
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use statistics
  
  use collection
  use convection
    
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
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ilvmin,ilvmax,nmdb
  
  ! Variable for a filename:
  character(LEN=SYS_STRLEN) :: sString
  character(LEN=SYS_STRLEN) :: sPRMFile, sTRIFile
  type(t_timer) :: rtimer
  ! make all the regions
  type(t_boundaryRegion), dimension(:),pointer :: p_rregion
  
  real(dp), dimension(:,:), pointer :: p_DbdyEdg
  
  real(DP) :: dx,dy
  integer :: iregions,iend,iindex
  integer, dimension(2) :: Isize

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
      
  !---------------------------------------------------------------------------------
      
  call tria_initStandardMeshFromRaw (&
      rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)
      
  ! get the number of boundary segments in
  ! the boundary component
  iend = boundary_igetNsegments(rproblem%rboundary,1)

  rproblem%isegments=iend

  allocate(p_rregion(iend))

  ! create the boundary regions
  do iregions=1,iend
  ! Idea: create the regions, check in which region the parameter is
  call boundary_createRegion(rproblem%rboundary, 1,&
                           iregions,p_rregion(iregions))
  end do
      

  nmdb=rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation%NMBD

  Isize = (/2,2*iend/)
  
  call storage_new ('cc_initParamTriang', 'edgesBdy',&
      Isize, ST_DOUBLE, &
      rproblem%h_DedgesAtBoundary, ST_NEWBLOCK_NOINIT)
  
  call storage_getbase_double2D(rproblem%h_DedgesAtBoundary,p_DbdyEdg)

  ! create the boundary regions
  do iregions=1,iend

    ! get the x,y coordinates for the current parameter value
    call boundary_getCoords(rproblem%rboundary, 1,&
                            p_rregion(iregions)%dminParam, dx, dy)
          
    iindex=iregions*2-1
    p_DbdyEdg(1,iindex) = dx
    p_DbdyEdg(2,iindex) = dy
                                    
    ! get the x,y coordinates for the current parameter value
    call boundary_getCoords(rproblem%rboundary, 1,&
                            p_rregion(iregions)%dmaxParam, dx, dy)

    p_DbdyEdg(1,iindex+1) = dx
    p_DbdyEdg(2,iindex+1) = dy
                        
  end do
                        
  deallocate(p_rregion)

  !-----------------------------------------------------------------------------------------------

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
    
  !    call storage_new ('tria_genEdgesAtBoundary2D', 'KMBD', &
  !          rproblem%RlevelInfo(1)%rtriangulation%NVBD, ST_INT, &
  !          rproblem%RlevelInfo(1)%rtriangulation%h_IedgesAtBoundary, ST_NEWBLOCK_NOINIT)
          
                          

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
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
