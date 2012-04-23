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
  use meshmodification
  
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
  integer :: i,ilvmin,ilvmax,iconvertToTriangleMesh
  real(DP) :: ddisturbMeshFactor
  
    ! Variable for a filename:
    character(LEN=SYS_STRLEN) :: sString
    character(LEN=SYS_STRLEN) :: sPRMFile, sTRIFile
    type(t_timer) :: rtimer

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
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sParametrisation',sPRMFile,bdequote=.true.)
                              
    call parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sMesh',sTRIFile,bdequote=.true.)
    
    call parlst_getvalue_int (rproblem%rparamList,'PARAMTRIANG',&
                              'iconvertToTriangleMesh',iconvertToTriangleMesh,0)
    
    call parlst_getvalue_double (rproblem%rparamList,'PARAMTRIANG',&
                                 'ddisturbMeshFactor',ddisturbMeshFactor,0.0_DP)
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rproblem%rboundary, sPRMFile)
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation, &
        sTRIFile, rproblem%rboundary)
        
    ! Probably we convert the mesh to a pure triangular mesh.
    if (iconvertToTriangleMesh .ne. 0) then
      call tria_rawGridToTri(rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation)
    end if

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
    
    if (ddisturbMeshFactor .gt. 0.0_DP) then
      ! Use the mesh disturbance to disturb the fine mesh -- and
      ! with that all coarser meshes as they are connected.
      call meshmod_disturbMesh (rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation,&
          ddisturbMeshFactor)
    end if
    
    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimeGridGeneration = &
      rproblem%rstatistics%dtimeGridGeneration + rtimer%delapsedReal

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
