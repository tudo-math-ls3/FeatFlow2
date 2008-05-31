!##############################################################################
!# ****************************************************************************
!# <name> ccinitdomain </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic initialisation routines for CC3D:
!# Initialisation of the domain control and the triangulations on all levels.
!#
!# 1.) cc_initDomain
!#     -> Read domain, read triangulation, refine the mesh
!#
!# 2.) cc_doneDomain
!#     -> Remove the meshes of all levels from the heap
!#
!# </purpose>
!##############################################################################

MODULE ccinitdomain

  USE fsystem
  USE storage
  USE linearsolver
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
  USE ccdomaincontrol
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initDomain (rproblem)
  
!<description>
  ! This routine initialises the domain control and triangulation of the
  ! domain. The corresponding .tri file is read from disc and
  ! the triangulation is refined as described by the NLMIN/NLMAX parameters
  ! from the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! Variable for a filename:  
  CHARACTER(LEN=SYS_STRLEN) :: sString
  CHARACTER(LEN=SYS_STRLEN) :: sTRIFile
  
  ! A pointer to the triangulation
  TYPE(t_triangulation), POINTER :: p_rtria
  
    ! First of all, initialise the domain control
    CALL ccdc_init(rproblem)

    ! Get a pointer to the coarse-most triangulation
    p_rtria => rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation

    ! Read in the .TRI filename from the parameter list
    CALL parlst_getvalue_string (rproblem%rparamList,'DOMAIN',&
                                   'sMesh',sString)
    READ (sString,*) sTRIFile

    ! And read the .TRI file
    CALL tria_readTriFile3D (p_rtria, sTRIFile)
    
    ! Call the domain control to correct the mesh (if neccessary)
    CALL ccdc_correctMesh(rproblem, p_rtria)

    ! Refine the mesh up to the minimum level
    DO i = 1, rproblem%NLMIN-1
      
      ! Refine it
      CALL tria_quickRefine2LevelOrdering(1, p_rtria)
      
      ! And call the mesh correction routine
      CALL ccdc_correctMesh(rproblem, p_rtria)
    
    END DO

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    CALL tria_initStandardMeshFromRaw (p_rtria)
    
    ! Now, refine to level up to nlmax.
    DO i=rproblem%NLMIN+1,rproblem%NLMAX
      
      ! Get a pointer to the triangulation
      p_rtria => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Refine the mesh
      CALL tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
                                      p_rtria)
      
      ! And call the mesh correction routine
      CALL ccdc_correctMesh(rproblem, p_rtria)
          
      ! And create all the information we need
      CALL tria_initStandardMeshFromRaw (p_rtria)
      
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

  SUBROUTINE cc_doneDomain (rproblem)
  
!<description>
  ! Releases the triangulation and domain control from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    ! Release the triangulation on all levels
    DO i = rproblem%NLMAX, rproblem%NLMIN,-1
      CALL tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! Release the domain control
    CALL ccdc_done(rproblem)
    
  END SUBROUTINE

END MODULE
