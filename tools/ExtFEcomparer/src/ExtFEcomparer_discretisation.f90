module ExtFEcomparer_discretisation

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use fparser
  use paramlist

  use triangulation
  use meshgeneration
  use basicgeometry

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

!<subroutine>

  subroutine ExtFEcomparer_initDynamicLevelInfo (rdynamicLevelInfo)

!<description>
  ! Initialises a dynamic level information structure with basic information.
!</description>

!<output>
  ! A dynamic level information structure to be initialised.
  type(t_dynamicLevelInfo), intent(out), target :: rdynamicLevelInfo
!</output>

!</subroutine>

    ! Initialise the BC/FBC structures.
    call bcasm_initDiscreteBC(rdynamicLevelInfo%rdiscreteBC)
    call bcasm_initDiscreteFBC(rdynamicLevelInfo%rdiscreteFBC)

    rdynamicLevelInfo%bhasNeumannBoundary = .false.
    rdynamicLevelInfo%hedgesDirichletBC = ST_NOHANDLE
    rdynamicLevelInfo%nedgesDirichletBC = 0

  end subroutine

  ! ***************************************************************************


 subroutine reinitDiscretisation_2D(rproblem)

!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>


!</subroutine>

  ! local variables
  integer :: i,j,ielementType,icubA,icubB,icubM
  integer :: iElementTypeStabil
  character(LEN=SYS_NAMELEN) :: sstr

    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary

    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation

    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Which discretisation is to use?
    ielementType = rproblem%ielemtype

    ! Now set up discrezisation structures on all levels:

    do i=rproblem%NLMIN,rproblem%NLMAX

      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation

      ! -----------------------------------------------------------------------
      ! Initialise discretisation structures for the spatial discretisation
      ! -----------------------------------------------------------------------

      ! Initialise the block discretisation according to the element specifier.
      call ExtFEcomparer_getDiscretisation_2D(ielementType,p_rdiscretisation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)

      ! Initialise the dynamic level information structure with basic information.
      call ExtFEcomparer_initDynamicLevelInfo (rproblem%RlevelInfo(i)%rdynamicInfo)

    end do

  end subroutine

  ! ***************************************************************************


  subroutine ExtFEcomparer_getDiscretisation_2D (ielementType,rdiscretisation,&
      rtriangulation,rboundary,rsourceDiscretisation)

!<description>
  ! Initialises a discretisation structure according to an element combination
  ! identifier.
!</description>

!<input>
  ! Element combination identifier to use. The following identifiers are supported:
  !  0 = Q1~(E031) / Q1~(E031) / Q0
  !  1 = Q1~(E030) / Q1~(E030) / Q0
  !  2 = Q1~(EM31) / Q1~(EM31) / Q0
  !  3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
  !  4 = Q2 (E013) / Q2 (E013) / QP1
  !  5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  !  6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  !  7 = Q1~(EM30) / Q1~(EM30) / Q0 new interface implementation
  !  8 = Q1~(EB30) / Q1~(EB30) / Q0
  ! 10 = Q2~(EB50) / Q2~(EB50) / QP1
  ! 11 = Q2~(EM50) / Q2~(EM50) / QP1
  ! 12 = Q2 (E013) / Q2 (E013) / QP1 (nonparametric)
  ! 20 = Q1        / Q1        / Q1
  ! 30 = P1~(E020) / P1~(E020) / P0
  integer, intent(in) :: ielementType

  ! Boundary structure for the discretisation
  type(t_boundary), intent(in), target, optional :: rboundary

  ! Triangulation structure for the discretisation
  type(t_triangulation), intent(in), target, optional :: rtriangulation

  ! OPTIONAL: Source discretisation. If specified, rdiscretisation is created
  ! as copy of rsourceDiscretisation with changed element type.
  ! If not specified, rboundary and rtriangulation must be present to
  ! specify mesh and boundary.
  type(t_blockDiscretisation), intent(in), target, optional :: rsourceDiscretisation
!</input>

!<output>
  ! Block discretisation structure that defines the discretisation according
  ! to ielementType.
  type(t_blockDiscretisation), intent(out), target :: rdiscretisation
!</output>

!</subroutine>

  ! Number of equations in our problem. velocity+velocity+pressure = 3
  integer, parameter :: nequations = 3

  ! local variables
  integer(I32) :: ieltypeUV, ieltypeP

    ! Initialise the element type identifiers according to ielementType
    select case (ielementType)
    case (0)
      ieltypeUV = EL_E031
      ieltypeP = EL_Q0

    case (1)
      ieltypeUV = EL_E030
      ieltypeP = EL_Q0

    case (2)
      ieltypeUV = EL_EM31
      ieltypeP = EL_Q0

    case (3)
      ieltypeUV = EL_EM30
      ieltypeP = EL_Q0

    case (4)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1

    case (5)
      ieltypeUV = EL_EM30_UNPIVOTED
      ieltypeP = EL_Q0

    case (6)
      ieltypeUV = EL_EM30_UNSCALED
      ieltypeP = EL_Q0

    case (7)
      ieltypeUV = EL_EM30_NEW
      ieltypeP = EL_Q0

    case (8)
      ieltypeUV = EL_EB30
      ieltypeP = EL_Q0

    case (10)
      ieltypeUV = EL_EB50
      ieltypeP = EL_QP1

    case (11)
      ieltypeUV = EL_EM50
      ieltypeP = EL_QP1

    case (12)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NP

    case (13)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NPD

    case (14)
      ieltypeUV = EL_Q1TB_2D
      ieltypeP = EL_QP1

    case (15)
      ieltypeUV = EL_Q1TBNP_2D
      ieltypeP = EL_QP1NPD

    case (20)
      ieltypeUV = EL_Q1
      ieltypeP = EL_Q1

    case (30)
      ieltypeUV = EL_P1T
      ieltypeP = EL_P0

    case (50)
      ieltypeUV = EL_EN50_2D
      ieltypeP  = EL_DCQP1_2D

    case (51)
      ieltypeUV = EL_EN51_2D
      ieltypeP  = EL_DCQP2_2D

    case (52)
      ieltypeUV = EL_QPW4P2_2D
      ieltypeP  = EL_QPW4DCP1_2D

    case default
      call output_line (&
          "Unknown discretisation: iElementType = "//sys_siL(ielementType,10), &
          OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D")
      call sys_halt()
    end select

    if (present(rsourceDiscretisation)) then

      ! Copy the source discretisation structure.
      call spdiscr_deriveBlockDiscr (rsourceDiscretisation,rdiscretisation)

      ! Replace the spatial discretisation structures inside by new ones
      ! with the correct element type.
      call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(1),  &
          ieltypeUV, rdiscretisation%RspatialDiscr(1))

      call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(2),  &
          ieltypeUV, rdiscretisation%RspatialDiscr(2))

      call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(3),  &
          ieltypeP, rdiscretisation%RspatialDiscr(3))

    else
      ! Ensure that rboundary/rtriangulation is present
      if ((.not. present(rboundary)) .or. &
          (.not. present(rtriangulation))) then
        call output_line (&
            "Boundary or triangulation not specified.", &
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D")
        call sys_halt()
      end if

      ! Check the mesh to ensure that the discretisation is compatible.
      select case(elem_igetShape(ieltypeUV))
      case (BGEOM_SHAPE_TRIA)
        ! Triangular element; ensure that we don"t have quads in the mesh
        if(rtriangulation%InelOfType(TRIA_NVEQUAD2D) .ne. 0) then
          call output_line (&
              "Discretisation does not support the current mesh!", &
              OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D")
          call sys_halt()
        end if

      case (BGEOM_SHAPE_QUAD)
        ! Quadrilateral element; ensure that we don"t have triangles in the mesh
        if(rtriangulation%InelOfType(TRIA_NVETRI2D) .ne. 0) then
          call output_line (&
              "Discretisation does not support the current mesh!", &
              OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D")
          call sys_halt()
        end if

      end select

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector.
      call spdiscr_initBlockDiscr (&
          rdiscretisation,nequations,rtriangulation,rboundary)

      ! rdiscretisation%RspatialDiscr is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
          rdiscretisation%RspatialDiscr(1), &
          ieltypeUV,rtriangulation, rboundary)

      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards.
      call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
          rdiscretisation%RspatialDiscr(2))

      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1),  &
          ieltypeP, rdiscretisation%RspatialDiscr(3))

    end if

  end subroutine


  subroutine ExtFEcomparer_doneDiscretisation (rproblem)

    !<description>
    ! Releases the discretisation from the heap.
    !</description>

    !<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
    !</inputoutput>

    !</subroutine>

    ! local variables
    integer :: i

    do i=rproblem%NLMAX,rproblem%NLMIN,-1

      ! Remove the block discretisation structure and all substructures.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation)
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisationStabil)
      call spdiscr_releaseDiscr(rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationStabil)

      ! Release dynamic level information
      call ExtFEcomparer_doneDynamicLevelInfo (rproblem%RlevelInfo(i)%rdynamicInfo)

    end do

end subroutine



  subroutine ExtFEcomparer_doneDynamicLevelInfo (rdynamicLevelInfo)

!<description>
  ! Releases a dynamic level information structure.
!</description>

!<inputoutput>
  ! A dynamic level information structure to be initialised.
  type(t_dynamicLevelInfo), intent(inout), target :: rdynamicLevelInfo
!</inputoutput>

!</subroutine>

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdynamicLevelInfo%rdiscreteBC)

    ! as well as the discrete version of the BC`s for fictitious boundaries
    call bcasm_releaseDiscreteFBC (rdynamicLevelInfo%rdiscreteFBC)

    ! Release the Dirichlet edges.
    if (rdynamicLevelInfo%hedgesDirichletBC .ne. ST_NOHANDLE) then
      call storage_free (rdynamicLevelInfo%hedgesDirichletBC)
    end if
    rdynamicLevelInfo%nedgesDirichletBC = 0

end subroutine

end module
