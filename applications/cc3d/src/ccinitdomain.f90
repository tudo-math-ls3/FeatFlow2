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

module ccinitdomain

  use fsystem
  use storage
  use linearsolver
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
  
  use collection
  use convection
    
  use ccbasic
  use ccdomaincontrol
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initDomain (rproblem)
  
!<description>
  ! This routine initialises the domain control and triangulation of the
  ! domain. The corresponding .tri file is read from disc and
  ! the triangulation is refined as described by the NLMIN/NLMAX parameters
  ! from the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
  ! Variable for a filename:
  character(LEN=SYS_STRLEN) :: sString
  character(LEN=SYS_STRLEN) :: sTRIFile
  
  ! A pointer to the triangulation
  type(t_triangulation), pointer :: p_rtria
  
    ! First of all, initialise the domain control
    call ccdc_init(rproblem)

    ! Get a pointer to the coarse-most triangulation
    p_rtria => rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation

    ! Read in the .TRI filename from the parameter list
    call parlst_getvalue_string (rproblem%rparamList,'DOMAIN',&
                                   'sMesh',sString)
    read (sString,*) sTRIFile

    ! And read the .TRI file
    call tria_readTriFile3D (p_rtria, sTRIFile)
    
    ! Call the domain control to correct the mesh (if neccessary)
    call ccdc_correctMesh(rproblem, p_rtria)

    ! Refine the mesh up to the minimum level
    do i = 1, rproblem%NLMIN-1
      
      ! Refine it
      call tria_quickRefine2LevelOrdering(1, p_rtria)
      
      ! And call the mesh correction routine
      call ccdc_correctMesh(rproblem, p_rtria)
    
    end do

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (p_rtria)
    
    ! Now, refine to level up to nlmax.
    do i=rproblem%NLMIN+1,rproblem%NLMAX
      
      ! Get a pointer to the triangulation
      p_rtria => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Refine the mesh
      call tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
                                      p_rtria)
      
      ! And call the mesh correction routine
      call ccdc_correctMesh(rproblem, p_rtria)
          
      ! And create all the information we need
      call tria_initStandardMeshFromRaw (p_rtria)
      
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

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneDomain (rproblem)
  
!<description>
  ! Releases the triangulation and domain control from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    ! Release the triangulation on all levels
    do i = rproblem%NLMAX, rproblem%NLMIN,-1
      call tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    ! Release the domain control
    call ccdc_done(rproblem)
    
  end subroutine

end module
