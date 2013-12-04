module ExtFEcomparer_parametrisation

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use fparser
  use paramlist

  use triangulation
  use meshgeneration

  use element
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation

  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  use vectorio
  use ExtFEcomparer_typedefs

  implicit none

contains

subroutine recreate_parametrisation_2D(rproblem)

  !<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  integer :: i,ilvmin,ilvmax

  ! Variable for a filename:
  character(LEN=SYS_STRLEN) :: sPRMFile, sTRIFile

  !get all parameters
  ilvmax = rproblem%NLMAX
  ilvmin = rproblem%NLMIN
  sPRMFile = rproblem%sPRMFile
  sTRIFile = rproblem%sTRIFile

    !---------------------------------------------------------------------!
    ! Now we have read in the parameters.
    ! We should create the parametrisation now. It is done in the way it
    ! is done in cc2d.
    !---------------------------------------------------------------------!

    ! Allocate memory for all the levels.
    allocate(rproblem%RlevelInfo(1:ilvmax))

    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rproblem%rboundary, sPRMFile)

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

    ! Compress the level hierarchy.
    ! Share the vertex coordinates of all levels, so the coarse grid coordinates
    ! are "contained" in the fine grid coordinates. The effect is:
    ! 1.) Save some memory
    ! 2.) Every change in the fine grid coordinates also affects the coarse
    !     grid coordinates and vice versa.
    do i=rproblem%NLMAX-1,rproblem%NLMIN,-1
      call tria_compress2LevelOrdHierarchy (rproblem%RlevelInfo(i+1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation)
    end do

end subroutine


subroutine ExtFEcomparer_doneParamTriang (rproblem)

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
