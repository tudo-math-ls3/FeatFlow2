module ExtFEcomparer_paramtriang

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use paramlist

  use triangulation
  use meshgeneration

  use element
  use cubature
  use spatialdiscretisation

  use ExtFEcomparer_typedefs

  implicit none

  private
  public :: ExtFEcomparer_recreate_ParamTriang
  public :: ExtFEcomparer_doneParamTriang

contains


! This is the function to be called from outside,
! it will branch its way to the right routine
subroutine ExtFEcomparer_recreate_ParamTriang(rproblem)
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

    select case (rproblem%iDimension)

    case (ExtFE_NDIM1)
        call ExtFEcomparer_recreate_ParamTriang_1D(rproblem)
    case (ExtFE_NDIM2)
        call ExtFEcomparer_recreate_ParamTriang_2D(rproblem)
    case default
        call output_line("At the moment the dimension of your problem is not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_recreate_parametrisation: dimension")
        call sys_halt()
    end select

end subroutine

subroutine ExtFEcomparer_recreate_ParamTriang_1D(rproblem)
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>
   integer :: i, ilvlmax

    ilvlmax = rproblem%NLMAX

  ! Read coarse mesh from TRI-file
  call tria_readTriFile1D(rproblem%rtriangulation,&
                rproblem%sTRIFile)


  ! Refine it to the maximum level
  call tria_quickRefine2LevelOrdering(ilvlmax-1,&
           rproblem%rtriangulation)

! Init a standard mesh
  call tria_initStandardMeshFromRaw(rproblem%rtriangulation)

end subroutine


subroutine ExtFEcomparer_recreate_ParamTriang_2D(rproblem)

  !<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  integer :: ilvmax

    ilvmax = rproblem%NLMAX

    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rproblem%rboundary, rproblem%sPRMFile)

   ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%rtriangulation, &
        rproblem%sTRIFile, rproblem%rboundary)

    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(ilvmax-1,&
        rproblem%rtriangulation,rproblem%rboundary)

      call tria_initStandardMeshFromRaw(rproblem%rtriangulation,&
          rproblem%rboundary)

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

    ! Release the triangulation on all levels
      call tria_done (rproblem%rtriangulation)


    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)

end subroutine


end module
