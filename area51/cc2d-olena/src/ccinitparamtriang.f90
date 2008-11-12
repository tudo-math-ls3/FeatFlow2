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

MODULE ccinitparamtriang

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  
  USE collection
  USE convection
    
  USE ccbasic
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initParamTriang (rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the NLMIN/NLMAX parameters
  ! from the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,ilvmin,ilvmax
  
    ! Variable for a filename:  
    CHARACTER(LEN=SYS_STRLEN) :: sString
    CHARACTER(LEN=SYS_STRLEN) :: sPRMFile, sTRIFile

    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)
    
    ! Get the .prm and the .tri file from the parameter list.
    ! note that parlst_getvalue_string returns us exactly what stands
    ! in the parameter file, so we have to apply READ to get rid of
    ! probable ''!
    CALL parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sParametrisation',sString)
    READ (sString,*) sPRMFile
                              
    CALL parlst_getvalue_string (rproblem%rparamList,'PARAMTRIANG',&
                                 'sMesh',sString)
    READ (sString,*) sTRIFile
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    CALL boundary_read_prm(rproblem%rboundary, sPrmFile)
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation, &
        sTRIFile, rproblem%rboundary)

    ! Refine the mesh up to the minimum level
    CALL tria_quickRefine2LevelOrdering(rproblem%NLMIN-1,&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    CALL tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)
    
    ! Now, refine to level up to nlmax.
    DO i=rproblem%NLMIN+1,rproblem%NLMAX
      CALL tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      CALL tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
    END DO
    
    ! Compress the level hierarchy.
    ! Share the vertex coordinates of all levels, so the coarse grid coordinates
    ! are 'contained' in the fine grid coordinates. The effect is:
    ! 1.) Save some memory
    ! 2.) Every change in the fine grid coordinates also affects the coarse
    !     grid coordinates and vice versa.
    DO i=rproblem%NLMAX-1,rproblem%NLMIN,-1
      CALL tria_compress2LevelOrdHierarchy (rproblem%RlevelInfo(i+1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    ! Release the triangulation on all levels
    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      CALL tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%rboundary)
    
  END SUBROUTINE

END MODULE
