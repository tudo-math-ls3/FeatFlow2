!<description>
! Purpose of all these functions is:
! We come here with the problem structure
! and can then recrate the parametrisation and
! the discretisation
!</description>

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
        call ExtFE_recreate_ParamTriang_1D(rproblem)
    case (ExtFE_NDIM2)
        call ExtFE_recreate_ParamTriang_2D(rproblem)
    case(ExtFE_NDIM3)
        call ExtFE_recreate_ParamTriang_3D(rproblem)
    case default
        call output_line("At the moment the dimension of your problem is not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_recreate_parametrisation: dimension")
        call sys_halt()
    end select

end subroutine

! In 1D, we have a TriFile and refine it.
! After that, we still have to init a Standard mesh from that.
! With this, we are done and have recreated the paramTriang for a
! 1D-mesh
! 1D-specific: The call to tria_readTriFile1D

subroutine ExtFE_recreate_ParamTriang_1D(rproblem)
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

! In 2D, we have a TriFile and a PRM-File
! With the prm we create the boundary-structure.
! Then we read in a TRI-File that we can refine as much
! as we have to.
! After that, we still have to init a Standard mesh from that.
! With this, we are done and have recreated the paramTriang for a
! 2D-mesh
! 2D-specific: The call to tria_readTriFile2D

subroutine ExtFE_recreate_ParamTriang_2D(rproblem)

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


! This routine just holds space for the 3D-specific routines
subroutine ExtFE_recreate_ParamTriang_3D(rproblem)
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

        call output_line("At the moment 3D is not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_recreate_parametrisation: dimension")
        call sys_halt()
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
